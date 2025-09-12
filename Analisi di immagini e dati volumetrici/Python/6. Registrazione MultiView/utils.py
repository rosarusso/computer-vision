import numpy as np
import open3d as o3d
from scipy.spatial.distance import cdist

def closestp_x84(data, model):
    # Compute distances between all points
    distances = cdist(data, model, metric='euclidean')
    
    # Find closest points
    min_indices = np.argmin(distances, axis=1)
    modelcp = model[min_indices]
    mindist = distances[np.arange(len(data)), min_indices]
    
    # Apply X84 strategy to discard outliers
    location = np.median(mindist)
    scale = 5.2 * np.median(np.abs(mindist - location))
    
    # Set outliers to NaN
    outlier_indices = np.abs(mindist - location) > scale
    modelcp[outlier_indices] = np.nan
    
    # Compute average distance for inliers
    inlier_indices = np.abs(mindist - location) <= scale
    res = np.mean(mindist[inlier_indices]) if np.any(inlier_indices) else 0
    
    return res, modelcp


def absolute(X, Y):
    # Discard NaN entries in X and corresponding entries in Y
    valid_indices = ~np.isnan(X).any(axis=1)
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
    
    # Create transformation matrix
    G = np.eye(4)
    G[:3, :3] = R
    G[:3, 3] = t
    
    return G


def eul(a):
    phi, theta, psi = a[2], a[1], a[0]
    
    Rz = np.array([
        [np.cos(phi), -np.sin(phi), 0],
        [np.sin(phi), np.cos(phi), 0],
        [0, 0, 1]
    ])
    
    Ry = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])
    
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(psi), -np.sin(psi)],
        [0, np.sin(psi), np.cos(psi)]
    ])
    
    R = Rz @ Ry @ Rx
    return R


def icp84(model, data):
    soglia = 1e-8
    G = np.eye(4)
    ris = np.inf
    risprev = 0
    i = 0
    
    while abs(ris - risprev) > soglia and i < 200:
        i += 1
        risprev = ris
        G = G[:3, :]  # Truncate to 3x4
        
        # Apply transformation
        dataHom = np.hstack([data, np.ones((data.shape[0], 1))])
        dataReg = (G @ dataHom.T).T
        
        # Find closest points
        distances = cdist(dataReg, model)
        min_indices = np.argmin(distances, axis=1)
        closest = model[min_indices]
        mindist = distances[np.arange(len(dataReg)), min_indices]
        
        # X84 outlier rejection
        MAD = np.median(np.abs(mindist - np.median(mindist)))
        out = 5.2 * MAD
        outlier_indices = np.abs(mindist - np.median(mindist)) > out
        closest[outlier_indices] = np.nan
        
        inlier_indices = np.abs(mindist - np.median(mindist)) <= out
        ris = np.mean(mindist[inlier_indices]) if np.any(inlier_indices) else 0
        
        # Discard NaNs
        valid_indices = ~np.isnan(closest).any(axis=1)
        closest = closest[valid_indices]
        dataReg = dataReg[valid_indices]
        
        # Compute centroids
        centroideX = np.mean(closest, axis=0)
        centroideY = np.mean(dataReg, axis=0)
        
        # Centered coordinates
        Xi = (closest - centroideX).T
        Yi = (dataReg - centroideY).T
        
        # SVD for rotation
        U, _, Vt = np.linalg.svd(Yi @ Xi.T)
        I = np.eye(3)
        I[2, 2] = np.linalg.det(Vt.T @ U.T)
        R = Vt.T @ I @ U.T
        t = centroideX - R @ centroideY
        
        # New transformation
        Gnew = np.eye(4)
        Gnew[:3, :3] = R
        Gnew[:3, 3] = t
        
        G = np.vstack([G, [0, 0, 0, 1]])
        G = Gnew @ G
    
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
        
        # Apply transformation
        dataHom = np.hstack([data, np.ones((data.shape[0], 1))])
        dataREG = (G @ dataHom.T).T[:, :3]
        
        # Compute closest points
        distances = cdist(dataREG, model)
        min_indices = np.argmin(distances, axis=1)
        modelCP = model[min_indices]
        mindist = distances[np.arange(len(dataREG)), min_indices]
        
        location = np.median(mindist)
        
        if scale is None:
            scale = 5.2 * np.median(np.abs(mindist - location))
        
        # Discard outliers
        outlier_indices = np.abs(mindist - location) > scale
        modelCP[outlier_indices] = np.nan
        
        inlier_indices = np.abs(mindist - location) <= scale
        res = np.mean(mindist[inlier_indices]) if np.any(inlier_indices) else 0
        
        # Compute incremental transformation
        GI = absolute(modelCP, dataREG)
        G = GI @ G
    
    print(f'Iterations: {i}')
    return G


def ieul(R):
    phi = np.arctan2(R[1, 0], R[0, 0])
    theta = np.arcsin(-R[2, 0])
    psi = np.arctan2(R[2, 1], R[2, 2])
    
    return np.array([psi, theta, phi])


def readply(filename):
    mesh = o3d.io.read_triangle_mesh(filename)
    assert mesh.has_vertices() and mesh.has_triangles(), "PLY file must have vertices and faces"
    return mesh


def rigid(G, M):
    if G.shape[1] == 1:  # Vector input
        R = eul(G[:3])
        t = G[3:]
        G = np.hstack([R, t.reshape(-1, 1)])
    else:  # Matrix input
        G = G[:3, :]
    
    HM = np.hstack([M, np.ones((M.shape[0], 1))])
    D = (G @ HM.T).T
    return D


class TriangularMesh:
    def __init__(self, verts, faces):
        self.updatemesh(verts, faces)
    
    def updatemesh(self, verts, faces):
        assert verts.shape[1] == 3, "Vertices must be Nx3"
        assert faces.shape[1] == 3, "Faces must be Mx3"
        self.verts = verts
        self.faces = faces
        self.n_verts = verts.shape[0]
        self.n_faces = faces.shape[0]
        assert np.all(faces >= 0) and np.all(faces < self.n_verts), "Face indices out of bounds"
    
    def getverts(self):
        return self.verts
    
    def getfaces(self):
        return self.faces
    
    def getnverts(self):
        return self.n_verts
    
    def getnfaces(self):
        return self.n_faces


def writeply(verts, faces, filename, mode='binary'):
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(verts)
    mesh.triangles = o3d.utility.Vector3iVector(faces)
    o3d.io.write_triangle_mesh(filename, mesh, write_ascii=(mode == 'ascii'))

