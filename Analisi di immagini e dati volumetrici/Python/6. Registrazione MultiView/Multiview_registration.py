import numpy as np
from utils import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots


G = np.eye(4)
Regdata = []

for i in range(1, 24):
    filename = f'./Bunny/bunny-range-{i:03d}.ply'
    mesh = o3d.io.read_triangle_mesh(filename)
    verts = np.asarray(mesh.vertices)
    
    # Subsample data
    Subdata = verts[::40]
    
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

# Visualize registered data using Plotly
fig = go.Figure()

# Add each point cloud as a separate trace
colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

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
        name=f'View {idx + 1}'
    ))

fig.update_layout(
    title='Multi-view Registered Point Clouds',
    scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z',
        aspectmode='data'
    ),
    width=800,
    height=600
)

fig.show()
