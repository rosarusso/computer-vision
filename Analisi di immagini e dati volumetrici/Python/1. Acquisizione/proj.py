# ------------------------------------------------------------
# proj.py
# ------------------------------------------------------------
# Perspective projection from 3‑D points to image pixel coordinates.
# ------------------------------------------------------------
import numpy as np

def proj(P: np.ndarray, points_3d: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the perspective projection of 3‑D points.

    Parameters
    ----------
    P : (3, 4) array
        Camera matrix (intrinsic * [R | t]).
    points_3d : (N, 3) array
        3‑D points expressed in the camera coordinate system.

    Returns
    -------
    u, v : (N,) int arrays
        Pixel coordinates (rounded to the nearest integer).
    """
    # Append homogeneous coordinate = 1
    homogeneous = np.hstack([points_3d, np.ones((points_3d.shape[0], 1))])  # (N,4)
    # Multiply by the projection matrix (3x4) → (N,3)
    h2d = homogeneous @ P.T                     # (N,3)

    # Normalise by the third coordinate
    u = np.rint(h2d[:, 0] / h2d[:, 2]).astype(int)
    v = np.rint(h2d[:, 1] / h2d[:, 2]).astype(int)

    return u, v

