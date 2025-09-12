import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist

def mesh_laplacian(vertex, face):
    nvertex = vertex.shape[0]
    nface = face.shape[0]
    
    print(f'MESH_LAPLACIAN: Calc Laplacian matrix for {nvertex:5d} vertices...')
    
    # Initialize edge matrix
    edge = np.zeros((nvertex, nvertex))
    
    # Compute edge lengths for each face
    for i in range(nface):
        # Get the three vertices of the face
        v1 = vertex[face[i, 0], :]
        v2 = vertex[face[i, 1], :]
        v3 = vertex[face[i, 2], :]
        
        # Compute edge lengths
        diff1 = np.linalg.norm(v1 - v2)
        diff2 = np.linalg.norm(v2 - v3)
        diff3 = np.linalg.norm(v3 - v1)
        
        # Fill symmetric edge matrix
        edge[face[i, 0], face[i, 1]] = diff1
        edge[face[i, 1], face[i, 0]] = diff1
        edge[face[i, 1], face[i, 2]] = diff2
        edge[face[i, 2], face[i, 1]] = diff2
        edge[face[i, 2], face[i, 0]] = diff3
        edge[face[i, 0], face[i, 2]] = diff3
    
    # Compute Laplacian matrix
    lap = np.zeros((nvertex, nvertex))
    
    for i in range(nvertex):
        # Find neighbors
        k = np.where(edge[i, :] > 0)[0]
        ni = len(k)
        
        if ni > 0:
            hi = np.mean(edge[i, k])
            invhi = np.mean(1.0 / edge[i, k])
            
            lap[i, i] = -(4/hi) * invhi
            lap[i, k] = (4/(hi*ni)) * (1.0 / edge[i, k])
    
    # Convert to sparse matrices
    edge_sparse = csr_matrix(edge)
    lap_sparse = csr_matrix(lap)
    
    return lap_sparse, edge_sparse

