import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.spatial.distance import cdist

def apply_rigid_transform(G, points):
    homogeneous = np.hstack([points, np.ones((points.shape[0], 1))])
    transformed = (G @ homogeneous.T).T
    return transformed[:, :3]

# Load and subsample
model = np.loadtxt('a4000007.cnn')
perm_m = np.random.permutation(model.shape[0])[:model.shape[0]//3]
modelperm = model[perm_m]

data = np.loadtxt('a4000001.cnn')
outliers = data[1000:1501] + np.random.randn(501, 3) * 50 + 300
data = np.vstack([data, outliers])
perm_d = np.random.permutation(data.shape[0])[:data.shape[0]//3]
dataperm = data[perm_d]

# Create subplot figure
fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}]],
    subplot_titles=("Before ICP", "After ICP")
)

# Add before ICP data
fig.add_trace(go.Scatter3d(
    x=model[:, 0], y=model[:, 1], z=model[:, 2],
    mode='markers',
    marker=dict(size=1, color='blue'),
    name='Model (Before)',
    legendgroup='before',
    showlegend=True
), row=1, col=1)

fig.add_trace(go.Scatter3d(
    x=data[:, 0], y=data[:, 1], z=data[:, 2],
    mode='markers',
    marker=dict(size=1, color='red'),
    name='Data (Before)',
    legendgroup='before',
    showlegend=True
), row=1, col=1)

# ICP with X84
threshold = 1e-8
G = np.eye(4)
res = np.inf
resprev = 0
i = 0

while abs(res - resprev) > threshold and i < 200:
    i += 1
    resprev = res

    dataReg = apply_rigid_transform(G, dataperm)

    # Find closest points
    distances = cdist(dataReg, modelperm)
    min_indices = np.argmin(distances, axis=1)
    mindist = distances[np.arange(len(dataReg)), min_indices]
    closest = modelperm[min_indices].copy()

    # X84 outlier rejection
    MAD = np.median(np.abs(mindist - np.median(mindist)))
    out = 5.2 * MAD

    outlier_indices = np.where(np.abs(mindist - np.median(mindist)) > out)[0]
    closest[outlier_indices] = np.nan

    inlier_indices = np.where(np.abs(mindist - np.median(mindist)) <= out)[0]
    res = np.mean(mindist[inlier_indices]) if len(inlier_indices) > 0 else 0

    # Remove NaN correspondences
    valid = ~np.isnan(closest).any(axis=1)
    closest_clean = closest[valid]
    dataReg_clean = dataReg[valid]

    # Compute centroids
    centroid_model = np.mean(closest_clean, axis=0)
    centroid_data = np.mean(dataReg_clean, axis=0)

    # Centralize
    Xi = closest_clean - centroid_model
    Yi = dataReg_clean - centroid_data

    # SVD
    U, _, Vt = np.linalg.svd(Yi.T @ Xi)
    I = np.eye(3)
    I[2, 2] = np.linalg.det(Vt.T @ U.T)

    R = Vt.T @ I @ U.T
    t = centroid_model - R @ centroid_data

    Gnew = np.eye(4)
    Gnew[:3, :3] = R
    Gnew[:3, 3] = t

    G = Gnew @ G

# Apply final transformation
final_data = apply_rigid_transform(G, data)

# Add after ICP data
fig.add_trace(go.Scatter3d(
    x=modelperm[:, 0], y=modelperm[:, 1], z=modelperm[:, 2],
    mode='markers',
    marker=dict(size=1, color='blue'),
    name='Model (After)',
    legendgroup='after',
    showlegend=True
), row=1, col=2)

fig.add_trace(go.Scatter3d(
    x=final_data[:, 0], y=final_data[:, 1], z=final_data[:, 2],
    mode='markers',
    marker=dict(size=1, color='red'),
    name='Data (After)',
    legendgroup='after',
    showlegend=True
), row=1, col=2)

# Update layout
fig.update_layout(
    title_text="ICP Registration Results",
    scene=dict(
        xaxis=dict(title='X'),
        yaxis=dict(title='Y'),
        zaxis=dict(title='Z')
    ),
    scene2=dict(
        xaxis=dict(title='X'),
        yaxis=dict(title='Y'),
        zaxis=dict(title='Z')
    ),
    height=600,
    width=1200
)

fig.show()

