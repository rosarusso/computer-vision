import numpy as np
from utils import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import open3d as o3d

G = np.eye(4)
Regdata = []
Originaldata = []  # Store original point clouds

for i in range(1, 24):
    filename = f'./Bunny/bunny-range-{i:03d}.ply'
    mesh = o3d.io.read_triangle_mesh(filename)
    verts = np.asarray(mesh.vertices)
    
    # Subsample data
    Subdata = verts[::40]
    
    # Store original data for subplot
    Originaldata.append(Subdata)
    
    if i == 1:
        data = Subdata
        Regdata.append(verts)
    else:
        model = data
        data = Subdata
        Gnew = icp84(model, data)
        Gloop = G @ Gnew
        
        # Transform full data
        vertsHom = np.hstack([verts, np.ones((verts.shape[0], 1))])
        datareg = (Gloop @ vertsHom.T).T[:, :3]
        
        G = Gloop
        Regdata.append(datareg)

# Visualize both original and registered data using subplots
fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}]],
    subplot_titles=('Original Point Clouds', 'Registered Point Clouds')
)

# Color palette for visualization
colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

# Add original point clouds to the left subplot (column 1)
for idx, pts in enumerate(Originaldata):
    color = colors[idx % len(colors)]
    fig.add_trace(go.Scatter3d(
        x=pts[:, 0],
        y=pts[:, 1],
        z=pts[:, 2],
        mode='markers',
        marker=dict(
            size=2,
            color=color,
            opacity=0.6
        ),
        name=f'Original View {idx + 1}',
        showlegend=False  # Hide from legend to avoid clutter
    ), row=1, col=1)

# Add registered point clouds to the right subplot (column 2)
for idx, pts in enumerate(Regdata):
    color = colors[idx % len(colors)]
    fig.add_trace(go.Scatter3d(
        x=pts[:, 0],
        y=pts[:, 1],
        z=pts[:, 2],
        mode='markers',
        marker=dict(
            size=2,
            color=color,
            opacity=0.6
        ),
        name=f'Registered View {idx + 1}'
    ), row=1, col=2)

# Update layout
fig.update_layout(
    title='Point Cloud Registration Comparison',
    scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z',
        aspectmode='data'
    ),
    scene2=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z',
        aspectmode='data'
    ),
    width=1200,
    height=600
)

fig.show()

