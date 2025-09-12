import numpy as np
import matplotlib.pyplot as plt

def absolute(X, Y):
    # Discard NaN entries in X and correspondingly in Y
    valid_indices = ~np.isnan(X)
    X_clean = X[valid_indices].reshape(-1, 3)
    Y_clean = Y[valid_indices].reshape(-1, 3)

    dime = Y_clean.shape[0]

    # Compute centroids
    cm = np.sum(Y_clean, axis=0) / dime
    cd = np.sum(X_clean, axis=0) / dime

    # Subtract centroids
    Yb = rigid(np.hstack([np.zeros(3), -cm]), Y_clean)
    Xb = rigid(np.hstack([np.zeros(3), -cd]), X_clean)

    # Compute rotation using SVD
    K = Xb.T @ Yb
    U, _, Vt = np.linalg.svd(K)
    S = np.diag([1, 1, np.linalg.det(U @ Vt.T)])
    R = U @ S @ Vt

    # Compute translation
    t = cd - R @ cm

    # Build rigid transformation matrix
    G = np.block([[R, t.reshape(3, 1)],
                  [np.zeros((1, 3)), 1]])
    return G


def eul(a):
    psi, theta, phi = a  # yaw, pitch, roll

    Rz = np.array([[np.cos(phi), -np.sin(phi), 0],
                   [np.sin(phi),  np.cos(phi), 0],
                   [0,            0,           1]])

    Ry = np.array([[np.cos(theta),  0, np.sin(theta)],
                   [0,              1, 0],
                   [-np.sin(theta), 0, np.cos(theta)]])

    Rx = np.array([[1, 0,             0],
                   [0, np.cos(psi), -np.sin(psi)],
                   [0, np.sin(psi),  np.cos(psi)]])

    R = Rz @ Ry @ Rx
    return R


def ieul(R):
    phi = np.arctan2(R[1, 0], R[0, 0])
    theta = np.arcsin(-R[2, 0])
    psi = np.arctan2(R[2, 1], R[2, 2])

    a = np.array([psi, theta, phi])
    return a


def rigid(G, M):
    if G.ndim == 1 or G.shape[1] == 1:
        # G is a vector: [rx, ry, rz, tx, ty, tz]
        R = eul(G[:3])
        t = G[3:]
        G_matrix = np.hstack([R, t.reshape(3, 1)])
    else:
        # G is a matrix
        G_matrix = G[:3, :]

    # Homogeneous coordinates
    HM = np.vstack([M.T, np.ones(M.shape[0])])
    D = (G_matrix @ HM).T
    return D



def vectarrow(p0, p1):
    if len(p0) == 3 and len(p1) == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], 'b-')

        # Arrow head
        p = p1 - p0
        alpha = 0.1
        beta = 0.1
        eps = 1e-10

        hu = [p1[0] - alpha * (p[0] + beta * (p[1] + eps)),
              p1[0],
              p1[0] - alpha * (p[0] - beta * (p[1] + eps))]
        hv = [p1[1] - alpha * (p[1] - beta * (p[0] + eps)),
              p1[1],
              p1[1] - alpha * (p[1] + beta * (p[0] + eps))]
        hw = [p1[2] - alpha * p[2],
              p1[2],
              p1[2] - alpha * p[2]]

        ax.plot(hu, hv, hw, 'b-')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.grid(True)
        plt.show()
    elif len(p0) == 2 and len(p1) == 2:
        plt.figure()
        plt.plot([p0[0], p1[0]], [p0[1], p1[1]], 'b-')

        # Arrow head
        p = p1 - p0
        alpha = 0.1
        beta = 0.1
        eps = 1e-10

        hu = [p1[0] - alpha * (p[0] + beta * (p[1] + eps)),
              p1[0],
              p1[0] - alpha * (p[0] - beta * (p[1] + eps))]
        hv = [p1[1] - alpha * (p[1] - beta * (p[0] + eps)),
              p1[1],
              p1[1] - alpha * (p[1] + beta * (p[0] + eps))]

        plt.plot(hu, hv, 'b-')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid(True)
        plt.show()
    else:
        raise ValueError("p0 and p1 must have the same dimension (2D or 3D)")

