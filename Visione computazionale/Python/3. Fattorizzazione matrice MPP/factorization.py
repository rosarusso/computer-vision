import numpy as np
import scipy.linalg
from scipy.io import loadmat

def krt(P):
    """
    Factorize a 3x4 projection matrix P into K, R, t components
    where P = K * [R | t]
    """
    # Extract the first 3x3 matrix from P
    P_3x3 = P[:, :3]
    
    # QR decomposition of the inverse of P_3x3
    Q, U = scipy.linalg.qr(np.linalg.inv(P_3x3))
    
    # Ensure U(3,3) > 0 (MATLAB indexing, so Python index [2,2])
    if U[2, 2] < 0:
        U = -1 * U
    
    # Ensure U(1,1) > 0 (MATLAB indexing, so Python index [0,0])
    if U[0, 0] < 0:
        E = np.array([[-1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        U = E @ U  # Dot product
        Q = Q @ E
    
    # Ensure U(2,2) > 0 (MATLAB indexing, so Python index [1,1])
    if U[1, 1] < 0:
        E = np.array([[1, 0, 0],
                      [0, -1, 0],
                      [0, 0, 1]])
        U = E @ U
        Q = Q @ E
    
    # Ensure det(Q) > 0
    if np.linalg.det(Q) < 0:
        Q = -Q
    
    # Compute scale factor
    s = np.linalg.det(Q)
    
    # Compute R and t
    R = s * Q.T
    t = s * U @ P[:, 3]
    
    # Compute K
    K = np.linalg.inv(U / U[2, 2])  # U(3,3) in MATLAB is U[2,2] in Python
    
    return K, R, t

# MPP factorization
# Rosa Russo VR445639

# Camera calibration:
# - obtain internal parameters (intrinsic to the camera itself),
#   which allow a mapping between camera coordinates and pixel coordinates
#   in the image frame;
# - obtain extrinsic parameters, that define the location and orientation
#   of the camera wrt the world frame.

# Load camera calibration data
camera = loadmat('Calib_Results.mat')

# Intrinsic parameters matrix
K = camera['KK']

# Rotation matrix
R = camera['Rc_5']

# Translation vector
t = camera['Tc_5']

# Construct projection matrix P
P = K @ np.column_stack((R, t))

# Factorize P back into K, R, t components
K_fact, R_fact, t_fact = krt(P)

# Print results
print("Original K:")
print(K)
print("\nFactorized K:")
print(K_fact)
print("\nOriginal R:")
print(R)
print("\nFactorized R:")
print(R_fact)
print("\nOriginal t:")
print(t)
print("\nFactorized t:")
print(t_fact)

