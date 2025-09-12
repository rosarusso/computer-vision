import numpy as np
from scipy.spatial.distance import cdist
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def apply_rigid_transform(G, points):
    """
    Apply a rigid transformation G (4x4 matrix) to a set of 3D points.
    
    Parameters:
        G (np.ndarray): 4x4 transformation matrix
        points (np.ndarray): Nx3 array of 3D points
    
    Returns:
        np.ndarray: Transformed Nx3 points
    """
    # Convert points to homogeneous coordinates
    ones = np.ones((points.shape[0], 1))
    points_homo = np.hstack([points, ones])
    
    # Apply transformation
    transformed = (G @ points_homo.T).T
    
    # Return back to 3D coordinates
    return transformed[:, :3]


def absolute(X, Y):
    # Discard rows with NaN in X or Y
    valid_indices = ~np.isnan(X).any(axis=1) & ~np.isnan(Y).any(axis=1)
    X = X[valid_indices]
    Y = Y[valid_indices]

    dime = Y.shape[0]

    # Compute centroids
    cm = np.mean(Y, axis=0)
    cd = np.mean(X, axis=0)

    # Subtract centroids
    Yb = Y - cm
    Xb = X - cd

    # Compute rotation using SVD
    K = Xb.T @ Yb
    U, _, Vt = np.linalg.svd(K)
    S = np.diag([1, 1, np.linalg.det(U @ Vt)])
    R = U @ S @ Vt

    # Compute translation
    t = cd - R @ cm

    # Build transformation matrix G
    G = np.eye(4)
    G[:3, :3] = R
    G[:3, 3] = t

    return G



def closestp(data, model):
    modelcp = np.zeros_like(data)
    mindist = np.full(data.shape[0], np.inf)

    for i in range(data.shape[0]):
        distances = np.linalg.norm(model - data[i], axis=1)
        min_idx = np.argmin(distances)
        mindist[i] = distances[min_idx]
        modelcp[i] = model[min_idx]

    res = np.mean(mindist)
    return res, modelcp



def closestp_x84(data, model):
    distances = cdist(data, model)
    min_indices = np.argmin(distances, axis=1)
    mindist = distances[np.arange(len(data)), min_indices]
    modelcp = model[min_indices].copy()

    location = np.median(mindist)
    scale = 5.2 * np.median(np.abs(mindist - location))

    # Identify outliers
    outlier_indices = np.where(np.abs(mindist - location) > scale)[0]

    # Set outliers to NaN
    modelcp[outlier_indices] = np.nan

    # Compute residual for inliers
    inlier_indices = np.where(np.abs(mindist - location) <= scale)[0]
    res = np.mean(mindist[inlier_indices]) if len(inlier_indices) > 0 else 0

    return res, modelcp


def eul(angles):
    psi, theta, phi = angles  # yaw, pitch, roll

    Rz = np.array([
        [np.cos(psi), -np.sin(psi), 0],
        [np.sin(psi),  np.cos(psi), 0],
        [0,            0,          1]
    ])

    Ry = np.array([
        [np.cos(theta),  0, np.sin(theta)],
        [0,              1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

    Rx = np.array([
        [1, 0,             0],
        [0, np.cos(phi), -np.sin(phi)],
        [0, np.sin(phi),  np.cos(phi)]
    ])

    R = Rz @ Ry @ Rx
    return R


def icp(model, data, maxiter=200):
    eps = 1e-8
    G = np.eye(4)
    res = np.inf
    prevres = 0
    i = 0

    while abs(res - prevres) > eps and i < maxiter:
        i += 1
        prevres = res

        dataREG = apply_rigid_transform(G, data)
        res, modelCP = closestp(dataREG, model)

        GI = absolute(modelCP, dataREG)
        G = GI @ G

    print(f"Iterations: {i}")
    return G


def icp_x84(model, data, maxiter=200, scale=None):
    eps = 1e-8
    G = np.eye(4)
    res = np.inf
    prevres = 0
    i = 0

    while abs(res - prevres) > eps and i < maxiter:
        i += 1
        prevres = res

        dataREG = apply_rigid_transform(G, data)
        res, modelCP = closestp_x84(dataREG, model)

        if scale is None:
            location = np.median(res)
            scale = 5.2 * np.median(np.abs(res - location))

        GI = absolute(modelCP, dataREG)
        G = GI @ G

    print(f"Iterations: {i}")
    return G


def ieul(R):
    psi = np.arctan2(R[2,1], R[2,2])  # yaw
    theta = np.arcsin(-R[2,0])        # pitch
    phi = np.arctan2(R[1,0], R[0,0])  # roll
    return np.array([psi, theta, phi])


def rigid(G, M):
    if G.ndim == 1:
        # Assume G is [roll, pitch, yaw, tx, ty, tz]
        R = eul(G[:3])
        t = G[3:]
        transformation_matrix = np.eye(4)
        transformation_matrix[:3, :3] = R
        transformation_matrix[:3, 3] = t
        G = transformation_matrix

    # Convert to homogeneous coordinates
    HM = np.hstack([M, np.ones((M.shape[0], 1))])
    D = (G @ HM.T).T
    return D[:, :3]

