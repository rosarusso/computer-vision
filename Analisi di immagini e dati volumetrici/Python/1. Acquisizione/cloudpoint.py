# ------------------------------------------------------------
# cloudpoint.py
# ------------------------------------------------------------
# Back‑project a depth image to a 3‑D point cloud and plot it.
# ------------------------------------------------------------

import pathlib
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio   # works with .pgm, .png, .tif, ...

# ------------------------------------------------------------------
# 1.  Intrinsic camera matrix 
# ------------------------------------------------------------------
K = np.array([[575.8, 0.0, 319.5],
              [0.0, 575.8, 239.5],
              [0.0,   0.0,   1.0]])

fx, fy = K[0, 0], K[1, 1]
cx, cy = K[0, 2], K[1, 2]

# ------------------------------------------------------------------
# 2.  Load depth image (the .pgm supplied with the repo)
# ------------------------------------------------------------------
depth_path = pathlib.Path('000_000595-b_8409599_depth.pgm')
if not depth_path.is_file():
    raise FileNotFoundError(f'Depth image not found: {depth_path}')
depth = imageio.imread(depth_path).astype(np.float64)   # shape (H, W)

# ------------------------------------------------------------------
# 3.  Back‑projection (vectorised – no Python loops!)
# ------------------------------------------------------------------
# Create a grid of pixel coordinates (x = column index, y = row index)
H, W = depth.shape
u = np.arange(W)            # columns 0 .. W-1
v = np.arange(H)            # rows    0 .. H-1
uu, vv = np.meshgrid(u, v)  # shape (H, W)

# Convert only non‑zero depth pixels
valid = depth != 0
x = (uu[valid] - cx) * depth[valid] / fx
y = (vv[valid] - cy) * depth[valid] / fy
z = depth[valid]

point_cloud = np.stack([x, y, z], axis=1)   # (N, 3)

# ------------------------------------------------------------------
# 4.  Visualisation – two sub‑plots
# ------------------------------------------------------------------
plt.figure(figsize=(12, 6))

# 4.1 Depth image
plt.subplot(1, 2, 1)
plt.imshow(depth, cmap='gray')
plt.title('Depth image')
plt.axis('off')

# 4.2 3‑D point cloud
ax = plt.subplot(1, 2, 2, projection='3d')
ax.scatter(point_cloud[:, 0],
           point_cloud[:, 1],
           point_cloud[:, 2],
           c='r', s=0.5)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Point cloud')
ax.view_init(elev=30, azim=-60)   # optional, nicer view
plt.show()

# ------------------------------------------------------------------
# 5.  Print the intrinsic matrix 
# ------------------------------------------------------------------
print('Perspective projection matrix:')
print(K)

