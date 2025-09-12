# ------------------------------------------------------------
# ply_io.py
# ------------------------------------------------------------

import numpy as np
from pathlib import Path
from typing import Iterable, Tuple, Union

def export_mesh_to_ply(
    vertices: np.ndarray,
    faces: np.ndarray,
    vertex_color: np.ndarray,
    name: Union[str, Path]
) -> None:
    """
    Write an ASCII PLY file.

    Parameters
    ----------
    vertices : (N, 3) float array
        XYZ coordinates of each vertex.
    faces    : (M, 3) int array
        0‑based vertex indices for each triangle.
    vertex_color : (N, 3) or (N,) uint8/float array
        Per‑vertex colour. If values are in [0, 1] they are scaled to
        [0, 255] automatically.
    name : str or pathlib.Path
        Output filename *without* the ".ply" extension (the function adds it).
    """
    vertices = np.asarray(vertices, dtype=np.float64)
    faces    = np.asarray(faces,    dtype=np.int32)

    # ------------------------------------------------------------------
    # 1. Normalise colour to uint8 
    # ------------------------------------------------------------------
    vc = np.asarray(vertex_color, dtype=np.float64)
    if vc.max() <= 1.0:
        vc = (vc * 256).clip(0, 255)
    if vc.ndim == 1:
        vc = np.tile(vc[:, None], (1, 3))
    vc = vc.astype(np.uint8)

    # ------------------------------------------------------------------
    # 2. Build the header
    # ------------------------------------------------------------------
    out_path = Path(name).with_suffix('.ply')
    with out_path.open('w') as f:
        f.write('ply\n')
        f.write('format ascii 1.0\n')
        f.write(f'element vertex {vertices.shape[0]}\n')
        f.write('property float x\n')
        f.write('property float y\n')
        f.write('property float z\n')
        f.write('property uchar red\n')
        f.write('property uchar green\n')
        f.write('property uchar blue\n')
        f.write(f'element face {faces.shape[0]}\n')
        f.write('property list uchar int vertex_index\n')
        f.write('end_header\n')

        # ------------------------------------------------------------------
        # 3. Write vertices
        # ------------------------------------------------------------------
        for (x, y, z), (r, g, b) in zip(vertices, vc):
            f.write(f'{x:.6f} {y:.6f} {z:.6f} {r} {g} {b}\n')

        # ------------------------------------------------------------------
        # 4. Write faces 
        # ------------------------------------------------------------------
        for v0, v1, v2 in faces:
            f.write(f'3 {v0} {v1} {v2}\n')

