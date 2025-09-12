
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.sparse.linalg import eigsh
from scipy.sparse import diags
import trimesh

def mesh_laplacian(vertex, triangle):
    n = vertex.shape[0]
    W = diags([1], [0], shape=(n, n), format='csr')
    A = diags([1], [0], shape=(n, n), format='csr')
    return W, A

def load_mesh(filepath):
    mesh = trimesh.load(filepath)
    vertex = mesh.vertices
    triangle = mesh.faces
    return vertex, triangle

# Load and compute eigenvalues for all three shapes
vertex1, triangle1 = load_mesh('./models/Male_scale.ply')
W1, A1 = mesh_laplacian(vertex1, triangle1)
val1, vet1 = eigsh(W1, 200, A1, sigma=-1e-5, which='LM')
val1 = np.abs(val1)

vertex2, triangle2 = load_mesh('./models/Male_null.ply')
W2, A2 = mesh_laplacian(vertex2, triangle2)
val2, vet2 = eigsh(W2, 200, A2, sigma=-1e-5, which='LM')
val2 = np.abs(val2)

vertex3, triangle3 = load_mesh('./models/Male_isometric.ply')
W3, A3 = mesh_laplacian(vertex3, triangle3)
val3, vet3 = eigsh(W3, 200, A3, sigma=-1e-5, which='LM')
val3 = np.abs(val3)

# Create subplots: 2 rows, 3 columns
fig = make_subplots(
    rows=2, cols=3,
    specs=[[{'type': 'mesh3d'}, {'type': 'mesh3d'}, {'type': 'mesh3d'}],
           [{'type': 'scatter'}, {'type': 'scatter'}, {'type': 'scatter'}]],
    subplot_titles=("Scale Shape", "Null Shape", "Isometric Shape",
                    "Scale Spectra", "Null Spectra", "Isometric Spectra")
)

# Add 3D meshes in the first row
fig.add_trace(go.Mesh3d(x=vertex1[:, 0], y=vertex1[:, 1], z=vertex1[:, 2],
                        i=triangle1[:, 0], j=triangle1[:, 1], k=triangle1[:, 2],
                        opacity=0.8, color='lightblue', name='Scale Mesh'),
              row=1, col=1)

fig.add_trace(go.Mesh3d(x=vertex2[:, 0], y=vertex2[:, 1], z=vertex2[:, 2],
                        i=triangle2[:, 0], j=triangle2[:, 1], k=triangle2[:, 2],
                        opacity=0.8, color='lightgreen', name='Null Mesh'),
              row=1, col=2)

fig.add_trace(go.Mesh3d(x=vertex3[:, 0], y=vertex3[:, 1], z=vertex3[:, 2],
                        i=triangle3[:, 0], j=triangle3[:, 1], k=triangle3[:, 2],
                        opacity=0.8, color='lightpink', name='Isometric Mesh'),
              row=1, col=3)

# Add eigenvalue spectra in the second row
fig.add_trace(go.Scatter(y=val1, mode='markers', marker=dict(color='green', symbol='x'), name='Scale Spectra'),
              row=2, col=1)

fig.add_trace(go.Scatter(y=val2, mode='markers', marker=dict(color='blue', symbol='circle'), name='Null Spectra'),
              row=2, col=2)

fig.add_trace(go.Scatter(y=val3, mode='markers', marker=dict(color='red', symbol='cross'), name='Isometric Spectra'),
              row=2, col=3)

# Update layout
fig.update_layout(
    title_text="Meshes and Eigenvalue Spectra Comparison",
    showlegend=False,
    height=800
)

# Show the combined figure
fig.show()

