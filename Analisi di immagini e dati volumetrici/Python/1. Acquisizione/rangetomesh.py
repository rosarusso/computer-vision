# ------------------------------------------------------------
# rangetomesh.py
# ------------------------------------------------------------
# Transform a depth/range image (with a binary validity mask) into a
# triangle mesh, applying the X‑84 outlier filter.
# ------------------------------------------------------------
import numpy as np
from tqdm import tqdm   # progress bar

def _distance(p1, p2):
    """Euclidean distance between two 3‑D points (both 1‑D arrays)."""
    return np.linalg.norm(p1 - p2)


def rangetomesh(valid_mask: np.ndarray,
                X: np.ndarray,
                Y: np.ndarray,
                Z: np.ndarray,
                K_X84: float = 5.2) -> tuple[np.ndarray, np.ndarray]:
    """
    Parameters
    ----------
    valid_mask : (H, W) bool or 0/1 array
        True where a depth sample exists.
    X, Y, Z : (H, W) float arrays
        The back‑projected coordinates
    K_X84 : float, optional
        Multiplier used in the MAD‑based outlier rejection (default 5.2).

    Returns
    -------
    triangles : (M, 3) int  (zero‑based vertex indices)
    vertices  : (N, 3) float
    """
    H, W = valid_mask.shape

    # ------------------------------------------------------------------
    # 1. Build dense vertex list and an index map (v_index) that tells
    #    for each pixel which vertex number it corresponds to.
    # ------------------------------------------------------------------
    v_index = np.zeros((H, W), dtype=np.int32) - 1   # -1 means "no vertex"
    vertices = []

    for i in range(H):
        for j in range(W):
            if valid_mask[i, j]:
                v_index[i, j] = len(vertices)   # 0‑based
                vertices.append([X[i, j], Y[i, j], Z[i, j]])

    vertices = np.asarray(vertices, dtype=np.float64)   # shape (N,3)

    # ------------------------------------------------------------------
    # 2. Compute the X‑84 distance threshold (TX_X84)
    # ------------------------------------------------------------------
    # Gather distances between a pixel and its 8‑connected neighbours.
    neighbour_offsets = [(-1, -1), (-1, 0), (-1, 1),
                         (0, -1),           (0, 1),
                         (1, -1),  (1, 0),  (1, 1)]

    dists = []
    for i in range(1, H - 1):
        for j in range(1, W - 1):
            if not valid_mask[i, j]:
                continue
            p = np.array([X[i, j], Y[i, j], Z[i, j]])
            for di, dj in neighbour_offsets:
                ii, jj = i + di, j + dj
                if valid_mask[ii, jj]:
                    q = np.array([X[ii, jj], Y[ii, jj], Z[ii, jj]])
                    dists.append(np.linalg.norm(p - q))
    dists = np.asarray(dists)

    # Median Absolute Deviation based threshold
    med = np.median(dists)
    mad = K_X84 * np.median(np.abs(dists - med))
    TX_X84 = med + mad

    # ------------------------------------------------------------------
    # 3. Build triangles 
    # ------------------------------------------------------------------
    triangles = []
    # Use tqdm for a progress bar (optional, can be removed)
    for i in tqdm(range(3, H - 3), desc='mesh generation', unit='row'):
        for j in range(3, W - 3):
            if not valid_mask[i, j]:
                continue

            # ---- central vertex -------------------------------------------------
            ind_central = v_index[i, j]
            central = vertices[ind_central]

            # ---- east neighbour -------------------------------------------------
            # check immediate east or fallback to east‑2 
            est_exists = False
            if valid_mask[i, j - 1]:
                est = vertices[v_index[i, j - 1]]
                if _distance(central, est) <= TX_X84:
                    est_idx = v_index[i, j - 1]
                    est_exists = True
            elif valid_mask[i, j - 2]:
                est = vertices[v_index[i, j - 2]]
                if _distance(central, est) <= TX_X84:
                    est_idx = v_index[i, j - 2]
                    est_exists = True

            # ---- north‑east neighbour -------------------------------------------
            nd_exists = False
            # 1) direct north‑east
            if valid_mask[i - 1, j - 1]:
                nd = vertices[v_index[i - 1, j - 1]]
                if _distance(central, nd) <= TX_X84:
                    nd_idx = v_index[i - 1, j - 1]
                    nd_exists = True
                    bool_nd = True
            # 2) fallbacks, respecting the “bool_est” flag
            if not nd_exists:
                # if we have a east neighbour (bool_est == True) we try (i‑1,j‑2)
                if est_exists and valid_mask[i - 1, j - 2]:
                    nd = vertices[v_index[i - 1, j - 2]]
                    if _distance(central, nd) <= TX_X84:
                        nd_idx = v_index[i - 1, j - 2]
                        nd_exists = True
                        bool_nd = True
                # else try the other combinations …
            if not nd_exists:
                if valid_mask[i - 2, j - 2]:
                    nd = vertices[v_index[i - 2, j - 2]]
                    if _distance(central, nd) <= TX_X84:
                        nd_idx = v_index[i - 2, j - 2]
                        nd_exists = True
                        bool_nd = False
            if not nd_exists:
                # try (i‑1,j‑2)
                if valid_mask[i - 1, j - 2]:
                    nd = vertices[v_index[i - 1, j - 2]]
                    if _distance(central, nd) <= TX_X84:
                        nd_idx = v_index[i - 1, j - 2]
                        nd_exists = True
                        bool_nd = False

            # ---- north neighbour -------------------------------------------------
            n_exists = False
            if valid_mask[i - 1, j]:
                n = vertices[v_index[i - 1, j]]
                if _distance(central, n) <= TX_X84:
                    n_idx = v_index[i - 1, j]
                    n_exists = True
            elif valid_mask[i - 2, j - 1]:
                n = vertices[v_index[i - 2, j - 1]]
                if _distance(central, n) <= TX_X84:
                    n_idx = v_index[i - 2, j - 1]
                    n_exists = True
            elif valid_mask[i, j - 2]:
                n = vertices[v_index[i, j - 2]]
                if _distance(central, n) <= TX_X84:
                    n_idx = v_index[i, j - 2]
                    n_exists = True

            # ---- Build the triangles  ----------------------
            # 1) (central, east, north‑east)
            if est_exists and nd_exists:
                if (_distance(central, vertices[est_idx]) <= TX_X84 and
                    _distance(central, vertices[nd_idx])  <= TX_X84 and
                    _distance(vertices[est_idx], vertices[nd_idx]) <= TX_X84):
                    triangles.append([ind_central, est_idx, nd_idx])

            # 2) (central, north‑east, north)
            if nd_exists and n_exists:
                if (_distance(central, vertices[nd_idx]) <= TX_X84 and
                    _distance(central, vertices[n_idx])   <= TX_X84 and
                    _distance(vertices[nd_idx], vertices[n_idx]) <= TX_X84):
                    triangles.append([ind_central, n_idx, nd_idx])

            # 3) fallback: (central, east, north) when north‑east missing
            if not nd_exists and n_exists and est_exists:
                if (_distance(central, vertices[est_idx]) <= TX_X84 and
                    _distance(central, vertices[n_idx])   <= TX_X84 and
                    _distance(vertices[est_idx], vertices[n_idx]) <= TX_X84):
                    triangles.append([ind_central, est_idx, n_idx])

    triangles = np.asarray(triangles, dtype=np.int32)
    return triangles, vertices

