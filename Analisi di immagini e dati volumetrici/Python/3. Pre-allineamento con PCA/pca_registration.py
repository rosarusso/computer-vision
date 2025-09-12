# pca_registration.py
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from utils import *

# Load data
data = np.loadtxt('a4000001.cnn')

# Generate random rotation vector
rotv = np.random.rand(3)

# Compute rotation matrices
Rz = np.array([[np.cos(rotv[2]), -np.sin(rotv[2]), 0],
               [np.sin(rotv[2]),  np.cos(rotv[2]), 0],
               [0,                0,               1]])

Ry = np.array([[np.cos(rotv[1]), 0, np.sin(rotv[1])],
               [0,               1, 0],
               [-np.sin(rotv[1]), 0, np.cos(rotv[1])]])

Rx = np.array([[1, 0,              0],
               [0, np.cos(rotv[0]), -np.sin(rotv[0])],
               [0, np.sin(rotv[0]),  np.cos(rotv[0])]])

R = Rz @ Ry @ Rx

# Build transformation matrix
G = np.hstack([R, np.zeros((3, 1))])

# Apply transformation
datamod = np.vstack([data.T, np.ones(data.shape[0])])
datarot = (G @ datamod).T

# Add noise
datarot += 20 * np.random.rand(*datarot.shape)

# Compute centroids
C1 = np.mean(data, axis=0)
C2 = np.mean(datarot, axis=0)

# PCA on data
_, _, vt1 = np.linalg.svd(data - C1)
u1 = vt1.T
R1 = u1[:, :3]

# Normalize PCA axes
for i in range(3):
    R1[:, i] /= np.linalg.norm(R1[:, i])

# Ensure right-hand system
if np.dot(np.cross(R1[:, 0], R1[:, 1]), R1[:, 2]) < 0:
    R1[:, 2] *= -1

# Apply transformation to data
G1 = np.hstack([R1.T, -R1.T @ C1.reshape(3, 1)])
datacan = rigid(np.hstack([np.zeros(3), -C1]), data)

# PCA on rotated data
_, _, vt2 = np.linalg.svd(datarot - C2)
u2 = vt2.T
R2 = u2[:, :3]

# Normalize PCA axes
for i in range(3):
    R2[:, i] /= np.linalg.norm(R2[:, i])

# Ensure right-hand system
if np.dot(np.cross(R2[:, 0], R2[:, 1]), R2[:, 2]) < 0:
    R2[:, 2] *= -1

# Generate 4 possible orientations
R2_all = np.zeros((3, 3, 4))
R2_all[:, :, 0] = R2
R2_all[:, :, 1] = np.column_stack([R2[:, 0], -R2[:, 1], -R2[:, 2]])
R2_all[:, :, 2] = np.column_stack([-R2[:, 0], R2[:, 1], -R2[:, 2]])
R2_all[:, :, 3] = np.column_stack([-R2[:, 0], -R2[:, 1], R2[:, 2]])

# Evaluate each case
Tr = np.zeros(4)
GReftot = np.zeros((4, 4, 4))

for i in range(4):
    R2 = R2_all[:, :, i]
    G2 = np.hstack([R2.T, -R2.T @ C2.reshape(3, 1)])
    G2mod = np.vstack([G2, [0, 0, 0, 1]])
    G1mod = np.vstack([G1, [0, 0, 0, 1]])
    GRef = np.linalg.inv(G1mod) @ G2mod
    GReftot[:, :, i] = GRef
    Tr[i] = np.trace(GRef[:3, :3])

# Select best case based on trace
optimal = np.argmax(Tr)
print(f"Best case according to trace heuristic: Case {optimal + 1}")
print(f"Trace value: {Tr[optimal]}")
print(f"Euler angles (deg): {np.rad2deg(ieul(GReftot[:3, :3, optimal]))}")

# Reconstruct the optimal transformation
R2_opt = R2_all[:, :, optimal]
G2_opt = np.hstack([R2_opt.T, -R2_opt.T @ C2.reshape(3, 1)])
datarotmod = np.vstack([datarot.T, np.ones(datarot.shape[0])])
datarot_can = (G2_opt @ datarotmod).T
GRef_opt = GReftot[:, :, optimal]
datarot_ref = (GRef_opt @ datarotmod).T

# Create a subplot figure with 1 row and 3 columns
fig = make_subplots(
    rows=1, cols=3,
    specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}, {'type': 'scatter3d'}]],
    subplot_titles=[
        "Original vs Rotated", "Canonical Pose", "Aligned Data"
    ],
    horizontal_spacing=0.02
)

# Add original vs rotated data
fig.add_trace(go.Scatter3d(x=data[:, 0], y=data[:, 1], z=data[:, 2],
                           mode='markers', marker=dict(size=1, color='red'),
                           name='Original Data'), row=1, col=1)
fig.add_trace(go.Scatter3d(x=datarot[:, 0], y=datarot[:, 1], z=datarot[:, 2],
                           mode='markers', marker=dict(size=1, color='blue'),
                           name='Rotated Data'), row=1, col=1)
fig.add_trace(go.Scatter3d(x=[C1[0]], y=[C1[1]], z=[C1[2]],
                           mode='markers', marker=dict(size=5, color='green'),
                           name='Centroid Original'), row=1, col=1)
fig.add_trace(go.Scatter3d(x=[C2[0]], y=[C2[1]], z=[C2[2]],
                           mode='markers', marker=dict(size=5, color='lightgreen'),
                           name='Centroid Rotated'), row=1, col=1)

# Add canonical pose subplot
fig.add_trace(go.Scatter3d(x=datacan[:, 0], y=datacan[:, 1], z=datacan[:, 2],
                           mode='markers', marker=dict(size=1, color='red'),
                           name='Canonical Original'), row=1, col=2)
fig.add_trace(go.Scatter3d(x=datarot_can[:, 0], y=datarot_can[:, 1], z=datarot_can[:, 2],
                           mode='markers', marker=dict(size=1, color='blue'),
                           name='Canonical Rotated'), row=1, col=2)

# Add aligned data subplot
fig.add_trace(go.Scatter3d(x=data[:, 0], y=data[:, 1], z=data[:, 2],
                           mode='markers', marker=dict(size=1, color='red'),
                           name='Original Data'), row=1, col=3)
fig.add_trace(go.Scatter3d(x=datarot_ref[:, 0], y=datarot_ref[:, 1], z=datarot_ref[:, 2],
                           mode='markers', marker=dict(size=1, color='blue'),
                           name='Aligned Rotated Data'), row=1, col=3)

# Update layout for better visualization
fig.update_layout(
    title=f"PCA Registration - Optimal Case {optimal + 1}",
    scene=dict(aspectmode='data'),
    height=600,
    width=1200
)

# Show the figure
fig.show()

