import numpy as np
from utils import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Load point clouds
model = np.loadtxt('a4000007.cnn')
data = np.loadtxt('a4000001.cnn')

# Subsample
model = model[np.random.permutation(model.shape[0])[:model.shape[0] // 3]]
data = data[np.random.permutation(data.shape[0])[:data.shape[0] // 3]]

# Create subplot figure with 3D plots
fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}]],
    subplot_titles=['Pre ICP', 'Post ICP']
)

# Plot before ICP
fig.add_trace(
    go.Scatter3d(
        x=model[:, 0],
        y=model[:, 1], 
        z=model[:, 2],
        mode='markers',
        marker=dict(color='blue', size=2),
        name='Model (Pre ICP)'
    ),
    row=1, col=1
)

fig.add_trace(
    go.Scatter3d(
        x=data[:, 0],
        y=data[:, 1],
        z=data[:, 2], 
        mode='markers',
        marker=dict(color='red', size=2),
        name='Data (Pre ICP)'
    ),
    row=1, col=1
)

# Run ICP
G = icp(model, data, maxiter=200)

# Apply final transformation
final_data = rigid(G, data)

# Plot after ICP
fig.add_trace(
    go.Scatter3d(
        x=model[:, 0],
        y=model[:, 1],
        z=model[:, 2],
        mode='markers',
        marker=dict(color='blue', size=2),
        name='Model (Post ICP)'
    ),
    row=1, col=2
)

fig.add_trace(
    go.Scatter3d(
        x=final_data[:, 0],
        y=final_data[:, 1],
        z=final_data[:, 2],
        mode='markers',
        marker=dict(color='red', size=2),
        name='Data (Post ICP)'
    ),
    row=1, col=2
)

# Set camera position and orientation
camera = dict(
    eye=dict(x=1.5, y=1.5, z=1.5),
    center=dict(x=0, y=0, z=0),
    up=dict(x=0, y=0, z=1)
)

# Update layout with camera settings
fig.update_layout(
    title_text='ICP Point Cloud Alignment',
    scene=dict(
        camera=camera,
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    ),
    scene2=dict(
        camera=camera,
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    ),
    showlegend=True
)

# Show the figure
fig.show()

