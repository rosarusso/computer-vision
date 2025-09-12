import numpy as np


def rigid(G, M):
    if G.ndim == 1 or G.shape[1] == 1:
        # G is a vector: [rx, ry, rz, tx, ty, tz]
        R = eul(G[:3])
        t = G[3:]
    else:
        # G is a 4x4 matrix
        R = G[:3, :3]
        t = G[:3, 3]

    # Apply transformation
    D = (R @ M.T).T + t
    return D


def ieul(R):
    phi = np.arctan2(R[1, 0], R[0, 0])
    theta = np.arcsin(-R[2, 0])
    psi = np.arctan2(R[2, 1], R[2, 2])

    return np.array([psi, theta, phi])


def icp(model, data, maxiter=200):
    eps = 1e-8
    G = np.eye(4)
    res = np.inf
    prevres = 0
    i = 0

    while abs(res - prevres) > eps and i < maxiter:
        i += 1
        prevres = res

        # Apply current transformation
        dataREG = rigid(G, data)

        # Find closest points
        _, modelCP = closestp(dataREG, model)

        # Compute residuals
        res = np.mean(np.linalg.norm(modelCP - dataREG, axis=1))

        # Compute incremental transformation
        GI = absolute(modelCP, dataREG)
        G = GI @ G

    print(f"Iterations: {i}")
    return G


def eul(a):
    psi, theta, phi = a  # yaw, pitch, roll

    Rz = np.array([
        [np.cos(phi), -np.sin(phi), 0],
        [np.sin(phi),  np.cos(phi), 0],
        [0,            0,           1]
    ])

    Ry = np.array([
        [np.cos(theta),  0, np.sin(theta)],
        [0,              1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

    Rx = np.array([
        [1, 0,             0],
        [0, np.cos(psi), -np.sin(psi)],
        [0, np.sin(psi),  np.cos(psi)]
    ])

    R = Rz @ Ry @ Rx
    return R



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


def absolute(X, Y):
    # Discard rows with NaN in X and corresponding rows in Y
    valid_indices = ~np.isnan(X).any(axis=1)
    X = X[valid_indices]
    Y = Y[valid_indices]

    dime = Y.shape[0]

    # Compute centroids
    cm = np.mean(Y, axis=0)
    cd = np.mean(X, axis=0)

    # Subtract centroids
    Yb = rigid(np.hstack(([0, 0, 0], -cm)), Y)
    Xb = rigid(np.hstack(([0, 0, 0], -cd)), X)

    # Compute rotation using SVD
    K = Xb.T @ Yb
    U, _, Vt = np.linalg.svd(K)
    S = np.diag([1, 1, np.linalg.det(U @ Vt)])
    R = U @ S @ Vt

    # Compute translation
    t = cd - R @ cm

    # Build transformation matrix
    G = np.eye(4)
    G[:3, :3] = R
    G[:3, 3] = t

    return G

