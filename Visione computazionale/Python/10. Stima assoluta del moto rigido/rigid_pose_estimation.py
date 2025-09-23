
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.io as sio

def eul(a):
    # Compute a rotation matrix given euler angles (yaw,pitch,roll)
    
    psi, theta, phi = a[0], a[1], a[2]
    
    # Rotation around z-axis (yaw)
    Rz = np.array([
        [np.cos(phi), -np.sin(phi), 0],
        [np.sin(phi),  np.cos(phi), 0],
        [0,            0,           1]
    ])
    
    # Rotation around y-axis (pitch)
    Ry = np.array([
        [np.cos(theta),  0, np.sin(theta)],
        [0,              1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])
    
    # Rotation around x-axis (roll)
    Rx = np.array([
        [1, 0,           0],
        [0, np.cos(psi), -np.sin(psi)],
        [0, np.sin(psi),  np.cos(psi)]
    ])
    
    R = Rz @ Ry @ Rx
    return R

def ieul(R):
    # Compute Euler angles (yaw,pitch,roll) given a rotation matrix
    
    phi = np.arctan2(R[1, 0], R[0, 0])
    theta = np.arcsin(-R[2, 0])
    psi = np.arctan2(R[2, 1], R[2, 2])
    
    a = np.array([psi, theta, phi])
    return a

def rigid(G, M):
    """ Apply 3D rigid transformation G to a point set M
    
    The transformation can be either a 4x4 homogeneous matrix or a 
    vector whose first 3 components are the rotation angles
    around the x,y,z axis respectively and the last 3 correspond to 
    the translation vector.
    """
	
    if G.ndim == 1: # Vector
        R = eul(G[:3])
        t = G[3:6].reshape(3, 1)
        G_matrix = np.hstack([R, t])
    else: # Matrix
        G_matrix = G[:3, :]
    
    # Convert M to homogeneous coordinates
    HM = np.hstack([M, np.ones((M.shape[0], 1))]).T
    
    # Apply transformation
    D = (G_matrix @ HM).T
    return D

def absolute(X, Y):
    """
    Solve absolute orientation 
    i.e. compute rigid transformation between two 3D point sets X and Y,
    When G is applied to Y, it brings Y on X. Entries containing NaN in X
    are discarded together with the corresponding entries in Y (used for 
    discarding correspondences)
    Algorithm ref.: Kanatani
    """
    # Flatten arrays to 1D for processing
    X_flat = X.flatten()
    Y_flat = Y.flatten()
    
    # Discard NaN entries in X and correspondingly in Y
    valid_indices = ~np.isnan(X_flat)
    X_valid = X_flat[valid_indices]
    Y_valid = Y_flat[valid_indices]
    
    # Reshape to Nx3 matrices
    dime = len(X_valid) // 3
    X = X_valid.reshape(dime, 3)
    Y = Y_valid.reshape(dime, 3)
    
    # Compute centroids
    cm = np.mean(Y, axis=0)
    cd = np.mean(X, axis=0)
    
    # Subtract centroids
    Yb = rigid(np.array([0, 0, 0, -cm[0], -cm[1], -cm[2]]), Y)
    Xb = rigid(np.array([0, 0, 0, -cd[0], -cd[1], -cd[2]]), X)
    
    # Compute rotation using SVD
    K = Xb.T @ Yb
    U, D, Vt = np.linalg.svd(K)
    V = Vt.T
    
    # Ensure proper rotation (determinant = 1)
    S = np.diag([1, 1, np.linalg.det(U @ V.T)])
    R = U @ S @ V.T
    
    # Compute translation
    t = cd.reshape(3, 1) - R @ cm.reshape(3, 1)
    
    # Create rigid transformation matrix
    G = np.block([[R, t], [np.array([0, 0, 0, 1])]])
    return G

# Load data
data = sio.loadmat('Corr3D.mat')
model_i = data['model_i']
data_i = data['data_i']

# Create first figure
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(model_i[:, 0], model_i[:, 1], model_i[:, 2], c='b', marker='.')
ax1.scatter(data_i[:, 0], data_i[:, 1], data_i[:, 2], c='r', marker='.')
ax1.set_title('Original Point Sets')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax1.grid(True)
plt.show()

# Compute rigid transformation between two 3D point sets "model_i" and "data_i"
G_out = absolute(model_i, data_i)

# Apply 3D rigid transformation G to a point set "data_i"
data_out = rigid(G_out, data_i)

# Create second figure
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(model_i[:, 0], model_i[:, 1], model_i[:, 2], c='b', marker='.')
ax2.scatter(data_out[:, 0], data_out[:, 1], data_out[:, 2], c='r', marker='.')
ax2.set_title('After Rigid Transformation')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
ax2.grid(True)
plt.show()

# Print transformation matrix
print("Rigid transformation matrix G_out:")
print(G_out)

