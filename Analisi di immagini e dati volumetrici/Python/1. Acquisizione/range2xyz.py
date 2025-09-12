# ------------------------------------------------------------
# range2xyz.py
# ------------------------------------------------------------
import trimesh
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from tqdm import tqdm   # optional, for the progress bar inside rangetomesh
from rangetomesh import rangetomesh
from ply_io import export_mesh_to_ply
import plotly.graph_objects as go
import numpy as np
from plotly.subplots import make_subplots

# ------------------------------------------------------------
# 1.  Load depth image
# ------------------------------------------------------------
depth_path = pathlib.Path('000_000595-b_8409599_depth.pgm')
if not depth_path.is_file():
    raise FileNotFoundError(f'Cannot find depth image: {depth_path}')
depth = imageio.imread(depth_path).astype(np.float64)   # shape (H, W)

# ------------------------------------------------------------
# 2.  Camera intrinsics
# ------------------------------------------------------------
fx, fy = 575.8, 575.8
cx, cy = 319.5, 239.5

# ------------------------------------------------------------
# 3.  Back‑project to 3‑D
# ------------------------------------------------------------
H, W = depth.shape
u = np.arange(W)
v = np.arange(H)
uu, vv = np.meshgrid(u, v)   # (H, W)

X = (uu - cx) * depth / fx
Y = (vv - cy) * depth / fy
Z = depth

# ------------------------------------------------------------
# 4.  Validity mask
# ------------------------------------------------------------
valid = depth != 0   # shape (H, W) bool

# ------------------------------------------------------------
# 5.  Build the mesh
# ------------------------------------------------------------
K_X84 = 5.2
triangles, vertices = rangetomesh(valid, X, Y, Z, K_X84)

# ------------------------------------------------------------
# 6.  Visualise the point cloud and the mesh
# ------------------------------------------------------------

# 3D scatter point cloud
scatter = go.Scatter3d(
    x=X[valid], y=Y[valid], z=Z[valid],
    mode='markers',
    marker=dict(color='red', size=2, opacity=0.6),
    name='Point cloud'
)

# Triangulated mesh
# Extract triangle vertex coordinates for Mesh3d
i, j, k = triangles.T  # triangles array shape (n_triangles, 3)

mesh = go.Mesh3d(
    x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
    i=i, j=j, k=k,
    color='lightgray',
    opacity=0.9,
    flatshading=True,
    name='Triangulated mesh',
    showscale=False,
    lighting=dict(ambient=0.5, diffuse=0.7, roughness=0.1, specular=0.2),
    lightposition=dict(x=100, y=200, z=0)
)


fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scene'}, {'type': 'scene'}]],
    subplot_titles=['Point cloud', 'Triangulated mesh']
)

# Add scatter to first subplot scene1
fig.add_trace(scatter, row=1, col=1)

# Add mesh to second subplot scene2
fig.add_trace(mesh, row=1, col=2)

camera = dict(
    eye=dict(x=0, y=0.1, z=1)
)

# Set axis titles for each subplot
fig.update_layout(
    scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z', camera=camera),
    scene2=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z', camera=camera),
    width=1300, height=900
)

fig.show()

