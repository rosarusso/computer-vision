import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from utils import *

# Read mesh
vertex, triangle = ply_read('./models/armadillo.ply')

# Plot original mesh
fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'mesh3d'}, {'type': 'mesh3d'}]])
fig.add_trace(go.Mesh3d(x=vertex[:, 0], y=vertex[:, 1], z=vertex[:, 2], 
                        i=triangle[:, 0], j=triangle[:, 1], k=triangle[:, 2], 
                        intensity=np.ones(len(vertex)), colorscale='Viridis', 
                        showscale=True), row=1, col=1)

# Compute Laplacian and normals
W, A = mesh_laplacian(vertex, triangle)
Am = np.diag(A)
L = np.linalg.inv(Am) @ W
# L = np.linalg.inv(A) @ W
N = getNormals(vertex, triangle)
curv = get_mcurvature(vertex, triangle, N, L)

# Plot curvature
fig.add_trace(go.Mesh3d(x=vertex[:, 0], y=vertex[:, 1], z=vertex[:, 2], 
                        i=triangle[:, 0], j=triangle[:, 1], k=triangle[:, 2], 
                        intensity=curv, colorscale='Viridis', showscale=True), 
             row=1, col=2)

fig.update_layout(title_text="Original Mesh and Mean Curvature")
fig.show()

