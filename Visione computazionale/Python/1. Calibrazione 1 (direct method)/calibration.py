
import cv2
import numpy as np
import matplotlib.pyplot as plt

# Calibration object coordinates (z, x, y, 1) - homogeneous coordinates
Mi = np.array([
    [0,   0,    0,   1],
    [0,   0,    8.3, 1],
    [5.3, 0,    8.3, 1],
    [5.3, 0,    0,   1],
    [0,   14.8, 8.3, 1],
    [5.3, 14.8, 8.3, 1]
])

# Load image
img1 = cv2.imread('Image1.jpg')
img1_rgb = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)  # Convert BGR to RGB for matplotlib

# Display image and collect points
fig, ax = plt.subplots()
ax.imshow(img1_rgb)
ax.set_title("Click on the 6 points in order")
plt.axis('off')

# Collect clicked points
mi = np.ones((len(Mi), 3))
points = plt.ginput(6, timeout=0)  # Wait indefinitely for 6 clicks

for i, (clickX, clickY) in enumerate(points):
    mi[i, :] = [clickX, clickY, 1]
    ax.scatter(clickX, clickY, c='green', marker='+')
    ax.text(clickX, clickY, f'. {i+1}', color='green')

plt.show()

# Calibration: build matrix A
A = np.zeros((2 * len(mi), 12))

for i in range(len(Mi)):
    ax = np.array([
        [0,       -mi[i, 2],  mi[i, 1]],
        [mi[i, 2],  0,       -mi[i, 0]],
        [-mi[i, 1], mi[i, 0],  0]
    ])
    KRO = np.kron(ax, Mi[i, :])
    A[2*i:2*i+2, :] = KRO[:2, :]

# Singular Value Decomposition
_, _, V = np.linalg.svd(A)
vecP = V[-1, :]  # Last column of V

# Reshape to get the Perspective Projection Matrix P
P = vecP.reshape(3, 4)

# Reproject points
m_reproj = np.zeros((len(Mi), 3))

for i in range(len(Mi)):
    mcurrent = P @ Mi[i, :].T
    m_reproj[i, :] = mcurrent / mcurrent[2]

# Display reprojected vs clicked points
fig, ax = plt.subplots()
ax.imshow(img1_rgb)
for i in range(len(Mi)):
    ax.scatter(m_reproj[i, 0], m_reproj[i, 1], c='red')
    ax.scatter(mi[i, 0], mi[i, 1], c='green', marker='+')
    ax.text(m_reproj[i, 0], m_reproj[i, 1], f'. {i+1}', color='red')
plt.axis('off')
plt.show()

# Save the projection matrix
img_i = 1
np.savez(f'Calib_direct_{img_i}.npz', P=P)

