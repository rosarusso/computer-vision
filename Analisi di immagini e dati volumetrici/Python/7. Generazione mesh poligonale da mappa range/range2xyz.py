import numpy as np
import matplotlib.pyplot as plt
import cv2
from utils import *


# Load depth image
DIZ = cv2.imread('000_000595-b_8409599_depth.pgm', cv2.IMREAD_UNCHANGED)

plt.figure(1)
plt.imshow(DIZ)

# Camera parameters
fx = 575.8
cx = 319.5
fy = 575.8
cy = 239.5

nrow, ncol = DIZ.shape
cloud = np.zeros((nrow*ncol, 3))
X = np.zeros((nrow, ncol))
Y = np.zeros((nrow, ncol))
TRUE = np.zeros((nrow, ncol))
k = 0
DephTH = 3000

# Generate point cloud
for i in range(ncol):  # column
    for j in range(nrow):  # row
        X[j,i] = (i - cx) * float(DIZ[j,i]) / fx
        Y[j,i] = (j - cy) * float(DIZ[j,i]) / fy
        if DIZ[j,i] != 0 and DIZ[j,i] < DephTH:
            TRUE[j,i] = 1
            k += 1
            cloud[k-1,:] = [X[j,i], Y[j,i], DIZ[j,i]]
    print(i)

# Plot point cloud
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(cloud[:,0], cloud[:,1], cloud[:,2], c='r', marker='.')

K_X84 = 5.2
Triangle, Vertex = rangetomesh(TRUE, X, Y, DIZ.astype(float), nrow, ncol, K_X84)

# Plot mesh
fig = plt.figure(3)
ax = fig.add_subplot(111, projection='3d')
# Note: You'll need to implement proper triangle mesh plotting

nvertex = Vertex.shape[0]
exportMeshToPly(Vertex, Triangle, np.ones((nvertex,3)), 'out2')

