import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.io import loadmat

data = loadmat('Corr3D.mat')
model_i = data['model_i']
data_i = data['data_i']

# Centroids computation
centroideX = np.mean(model_i, axis=0)
centroideY = np.mean(data_i, axis=0)

# Centralized coordinates
Xi = (model_i - centroideX).T
Yi = (data_i - centroideY).T

# SVD
U, S, Vt = np.linalg.svd(Yi @ Xi.T)

# Ensure proper rotation (det = 1)
I = np.eye(3)
I[2, 2] = np.linalg.det(Vt.T @ U)

# Rotation and translation
R = Vt.T @ I @ U.T
t = centroideX.T - R @ centroideY.T

# Transformation matrix
G = np.hstack([R, t.reshape(3, 1)])

# Apply transformation
data_imod = np.hstack([data_i, np.ones((data_i.shape[0], 1))]).T
data_out = (G @ data_imod).T

camera = dict(
    eye=dict(x=1, y=1, z=1)
)
# Create subplots with 3D scenes
fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}]],
    subplot_titles=('Before Registration', 'After Registration')
)

# Before Registration plot
fig.add_trace(
    go.Scatter3d(
        x=model_i[:, 0],
        y=model_i[:, 1], 
        z=model_i[:, 2],
        mode='markers',
        marker=dict(color='blue', size=3),
        name='Model (Before)'
    ),
    row=1, col=1
)

fig.add_trace(
    go.Scatter3d(
        x=data_i[:, 0],
        y=data_i[:, 1],
        z=data_i[:, 2],
        mode='markers',
        marker=dict(color='red', size=3),
        name='Data (Before)'
    ),
    row=1, col=1
)

# After Registration plot
fig.add_trace(
    go.Scatter3d(
        x=model_i[:, 0],
        y=model_i[:, 1],
        z=model_i[:, 2],
        mode='markers',
        marker=dict(color='blue', size=3),
        name='Model (After)'
    ),
    row=1, col=2
)

fig.add_trace(
    go.Scatter3d(
        x=data_out[:, 0],
        y=data_out[:, 1],
        z=data_out[:, 2],
        mode='markers',
        marker=dict(color='red', size=3),
        name='Data (After)'
    ),
    row=1, col=2
)

# Update layout
fig.update_layout(
    title_text="3D Point Cloud Registration",
    scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z', camera=camera),
    scene2=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z', camera=camera),
    showlegend=True,
    height=600,
    width=1200
)

# Update 3D scene layouts
fig.update_scenes(
    xaxis_title='X',
    yaxis_title='Y', 
    zaxis_title='Z',
    aspectmode='data'
)

fig.show()

print("Rototranslation matrix:")
print(G)

