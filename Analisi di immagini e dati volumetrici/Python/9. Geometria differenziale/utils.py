import numpy as np
from scipy.sparse import csr_matrix
from plyfile import PlyData
import time

def get_mcurvature(V, T, N, L):
    # Compute mean curvature
    m_curv = 0.5 * np.sum(N * (L @ V), axis=1)
    return m_curv


def getNormals(V, T):
    # Compute normals for each triangle
    normals = np.cross(V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 0]])
    
    # Accumulate normals per vertex
    N = np.zeros_like(V)
    for i in range(3):
        np.add.at(N, T[:, i], normals)
    
    # Normalize normals
    N /= np.linalg.norm(N, axis=1, keepdims=True)
    return N


def mesh_laplacian(vertex, face):
    """
    Compute the Laplacian matrix and edge distances for an irregular triangular mesh.

    Parameters:
    - vertex: np.ndarray of shape (N, 3), where N is the number of vertices.
               Each row contains the (x, y, z) Cartesian coordinates of a vertex.
    - face: np.ndarray of shape (M, 3), where M is the number of faces (triangles).
            Each row contains indices into 'vertex' that define a triangle.

    Returns:
    - lap: scipy.sparse.csr_matrix, the Laplacian matrix (N x N)
    - edge: scipy.sparse.csr_matrix, the edge distance matrix (N x N)

    Based on:
    Oostendorp, Oosterom & Huiskamp (1989),
    Interpolation on a triangulated 3D surface.
    Journal of Computational Physics, 80: 331-343.
    """

    nvertex = vertex.shape[0]
    nface = face.shape[0]

    print(f'MESH_LAPLACIAN: Calc Laplacian matrix for {nvertex:5d} vertices...', end='')

    start_time = time.time()

    # Initialize edge matrix as dense
    edge = np.zeros((nvertex, nvertex))

    # Compute edge lengths and populate edge matrix
    for i in range(nface):
        v1, v2, v3 = face[i]

        # Compute differences between vertices
        diff1 = vertex[v1] - vertex[v2]
        diff2 = vertex[v2] - vertex[v3]
        diff3 = vertex[v3] - vertex[v1]

        # Compute norms (distances)
        norm1 = np.linalg.norm(diff1)
        norm2 = np.linalg.norm(diff2)
        norm3 = np.linalg.norm(diff3)

        # Fill in symmetric edge distances
        edge[v1, v2] = norm1
        edge[v2, v1] = norm1

        edge[v2, v3] = norm2
        edge[v3, v2] = norm2

        edge[v3, v1] = norm3
        edge[v1, v3] = norm3

    # Initialize Laplacian matrix as dense
    lap = np.zeros((nvertex, nvertex))

    # Compute Laplacian values
    for i in range(nvertex):
        # Find neighbors: non-zero entries in the i-th row of edge
        neighbors = np.nonzero(edge[i, :])[0]
        ni = len(neighbors)

        if ni > 0:
            hi = np.mean(edge[i, neighbors])
            inv_hi = np.mean(1.0 / edge[i, neighbors])

            lap[i, i] = -(4 / hi) * inv_hi  # Self Laplacian

            # Neighbors Laplacian
            lap[i, neighbors] = (4 / (hi * ni)) * (1.0 / edge[i, neighbors])

    areas = np.zeros(nvertex)
    for i in range(nface):
        v1, v2, v3 = face[i]
        # Triangle area using cross product
        tri_area = 0.5 * np.linalg.norm(np.cross(vertex[v2] - vertex[v1], vertex[v3] - vertex[v1]))
        areas[v1] += tri_area / 3.0
        areas[v2] += tri_area / 3.0
        areas[v3] += tri_area / 3.0

    # Return the Laplacian matrix and vertex areas
    edge_sparse = csr_matrix(edge)
    lap_sparse = csr_matrix(lap)
    elapsed_time = time.time() - start_time
    print(f'done ({elapsed_time:6.2f} sec).')
    return lap_sparse, areas


def mshlp_matrix(shape, opt=None):
    if opt is None:
        opt = {}
    opt.setdefault('hs', 2)
    opt.setdefault('rho', 3)
    opt.setdefault('htype', 'ddr')
    opt.setdefault('dtype', 'cotangent')
    
    # Placeholder for meshlp function (needs to be implemented or imported)
    II, JJ, SS, AA = meshlp(shape['TRIV'], shape['X'], shape['Y'], shape['Z'], opt)
    
    W = csr_matrix((SS, (II, JJ)), shape=(len(shape['X']), len(shape['X'])))
    A = np.array(AA)
    
    return W, A


def ply_read(path):
    plydata = PlyData.read(path)
    
    # Extract vertex data
    vertex = plydata['vertex']
    vertex_data = np.column_stack((vertex['x'], vertex['y'], vertex['z']))
    
    # Extract face data
    face_data = [face for face in plydata['face']['vertex_indices']]
    face_indices = np.array(face_data)
    
    return vertex_data, face_indices

