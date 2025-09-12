import numpy as np
from tqdm import tqdm

def exportMeshToPly(vertices, faces, vertex_color, name):
    # Handle vertex color scaling
    if np.max(vertex_color) <= 1.0:
        vertex_color = vertex_color * 256
    
    # Handle single channel colors
    if vertex_color.shape[1] == 1:
        vertex_color = np.tile(vertex_color, (1, 3))
    
    vertex_color = vertex_color.astype(np.uint8)
    
    # Write PLY file
    with open(f"{name}.ply", 'w') as fidply:
        # Write header
        fidply.write('ply\n')
        fidply.write('format ascii 1.0\n')
        fidply.write(f'element vertex {vertices.shape[0]}\n')
        fidply.write('property float x\n')
        fidply.write('property float y\n')
        fidply.write('property float z\n')
        fidply.write('property uchar red\n')
        fidply.write('property uchar green\n')
        fidply.write('property uchar blue\n')
        fidply.write(f'element face {faces.shape[0]}\n')
        fidply.write('property list uchar int vertex_index\n')
        fidply.write('end_header\n')
        
        # Write vertices
        for i in range(vertices.shape[0]):
            fidply.write(f'{vertices[i,0]} {vertices[i,1]} {vertices[i,2]} {vertex_color[i,0]} {vertex_color[i,1]} {vertex_color[i,2]}\n')
        
        # Write faces (adjust indices to 0-based)
        for i in range(faces.shape[0]):
            fidply.write(f'3 {faces[i,0]} {faces[i,1]} {faces[i,2]}\n')


def proj(P, c3d):
    """
    PROJ : compute perspective projection (from 3D to pixel coordinates)
    """
    # Convert 3D coordinates to homogeneous coordinates
    h3d = np.hstack([c3d, np.ones((c3d.shape[0], 1))]).T
    
    # Apply projection matrix
    h2d = P @ h3d
    
    # Convert back to 2D coordinates
    c2d = h2d[:2,:] / h2d[2,:]
    
    # Round to pixel coordinates
    u = np.round(c2d[0,:]).astype(int)
    v = np.round(c2d[1,:]).astype(int)
    
    return u, v


def rangetomesh(TRUE, XX, YY, ZZ, dim_I, dim_J, K_X84):
    """
    Script per generare una mesh da una mappa range
    Nota che TRUE prende valori 0 o 1
    applica la X84 per eliminare gli edge lunghi
    """
    
    # Create vertex index structure
    v_index = np.ones((dim_I, dim_J))
    
    Vertex = np.zeros((dim_I*dim_J, 3))
    ind = 0
    
    # Create vertices
    for i in range(dim_I):
        for j in range(dim_J):
            if TRUE[i,j] == 1:
                ind += 1
                v_index[i,j] = ind
                Vertex[ind-1,:] = [XX[i,j], YY[i,j], ZZ[i,j]]
    
    Vertex = Vertex[:ind,:]
    
    # Calculate distance threshold using X84 rule
    res = np.zeros((dim_I*dim_J)*8)
    ind = 0
    
    for i in range(2, dim_I):
        for j in range(2, dim_J):
            # Check 8 neighbors in north-east quadrant
            if TRUE[i,j] == 1:
                v1 = np.array([XX[i,j], YY[i,j], ZZ[i,j]])
                for h in range(3):
                    for k in range(3):
                        if TRUE[i-h,j-k] == 1:
                            ind += 1
                            v2 = np.array([XX[i-h,j-k], YY[i-h,j-k], ZZ[i-h,j-k]])
                            res[ind-1] = np.linalg.norm(v1-v2)
    
    res = res[:ind]
    
    # Apply X84 rule: |x_i-MED|<k*MAD
    MED_x84 = np.median(res)
    MAD_x84 = K_X84 * np.median(np.abs(res - MED_x84))
    TX_X84 = MED_x84 + MAD_x84
    
    # Create triangles
    Triangle = np.zeros((dim_I*dim_J*8, 3))
    ind = 0
    
    # Progress bar equivalent
    for i in tqdm(range(2, dim_I-3), desc='mesh generation'):
        for j in range(2, dim_J-3):
            # Check that the point exists
            if TRUE[i,j] == 1:
                ind_central = int(v_index[i,j])
                central = Vertex[ind_central-1,:]  # Adjust for 0-based indexing
                
                # Neighbor 1: East
                est1 = Vertex[int(v_index[i,j-1])-1,:]
                diff1 = np.linalg.norm(central-est1)
                est2 = Vertex[int(v_index[i,j-2])-1,:]
                diff2 = np.linalg.norm(central-est2)
                
                if TRUE[i,j-1] == 1 and diff1 <= TX_X84:
                    ind_est = int(v_index[i,j-1])
                    bool_est = 1
                    esist_est = 1
                elif TRUE[i,j-2] == 1 and diff2 <= TX_X84:
                    ind_est = int(v_index[i,j-2])
                    bool_est = 0
                    esist_est = 1
                else:
                    esist_est = 0
                    bool_est = 1
                
                # Neighbor 2: North-East
                nordest1 = Vertex[int(v_index[i-1,j-1])-1,:]
                diff1 = np.linalg.norm(central-nordest1)
                nordest2 = Vertex[int(v_index[i-1,j-2])-1,:]
                diff2 = np.linalg.norm(central-nordest2)
                nordest3 = Vertex[int(v_index[i-2,j-2])-1,:]
                diff3 = np.linalg.norm(central-nordest3)
                
                if TRUE[i-1,j-1] == 1 and diff1 <= TX_X84:
                    ind_nordest = int(v_index[i-1,j-1])
                    bool_nordest = 1
                    esist_nordest = 1
                elif bool_est == 0 and TRUE[i-1,j-2] == 1 and diff2 <= TX_X84:
                    ind_nordest = int(v_index[i-1,j-2])
                    bool_nordest = 1
                    esist_nordest = 1
                elif bool_est == 0 and TRUE[i-2,j-2] == 1 and diff3 <= TX_X84:
                    ind_nordest = int(v_index[i-2,j-2])
                    bool_nordest = 0
                    esist_nordest = 1
                elif bool_est == 1 and TRUE[i-2,j-2] == 1 and diff3 <= TX_X84:
                    ind_nordest = int(v_index[i-2,j-2])
                    bool_nordest = 0
                    esist_nordest = 1
                elif bool_est == 1 and TRUE[i-1,j-2] == 1 and diff2 <= TX_X84:
                    ind_nordest = int(v_index[i-1,j-2])
                    bool_nordest = 0
                    esist_nordest = 1
                else:
                    esist_nordest = 0
                    bool_nordest = 1
                
                # Neighbor 3: North
                nord1 = Vertex[int(v_index[i-1,j])-1,:]
                diff1 = np.linalg.norm(central-nord1)
                nord2 = Vertex[int(v_index[i-2,j-1])-1,:]
                diff2 = np.linalg.norm(central-nord2)
                nord3 = Vertex[int(v_index[i,j-2])-1,:]
                diff3 = np.linalg.norm(central-nord3)
                
                if TRUE[i-1,j] == 1 and diff1 <= TX_X84:
                    ind_nord = int(v_index[i-1,j])
                    esist_nord = 1
                elif TRUE[i-2,j-1] == 1 and diff2 <= TX_X84:
                    ind_nord = int(v_index[i-2,j-1])
                    esist_nord = 1
                elif TRUE[i,j-2] == 1 and diff3 <= TX_X84:
                    ind_nord = int(v_index[i,j-2])
                    esist_nord = 1
                else:
                    esist_nord = 0
                
                # Create triangles
                # First triangle
                if esist_est == 1 and esist_nordest == 1:
                    if (np.linalg.norm(Vertex[ind_central-1,:]-Vertex[ind_est-1,:]) <= TX_X84 and 
                        np.linalg.norm(Vertex[ind_central-1,:]-Vertex[ind_nordest-1,:]) <= TX_X84 and 
                        np.linalg.norm(Vertex[ind_est-1,:]-Vertex[ind_nordest-1,:]) <= TX_X84):
                        ind += 1
                        Triangle[ind-1,:] = [ind_central-1, ind_est-1, ind_nordest-1]  # Convert to 0-based indexing
                
                # Second triangle
                if esist_nordest == 1 and esist_nord == 1:
                    if (np.linalg.norm(Vertex[ind_central-1,:]-Vertex[ind_nordest-1,:]) <= TX_X84 and 
                        np.linalg.norm(Vertex[ind_central-1,:]-Vertex[ind_nord-1,:]) <= TX_X84 and 
                        np.linalg.norm(Vertex[ind_nord-1,:]-Vertex[ind_nordest-1,:]) <= TX_X84):
                        ind += 1
                        Triangle[ind-1,:] = [ind_central-1, ind_nord-1, ind_nordest-1]  # Convert to 0-based indexing
                
                # If north-east doesn't exist, generate triangle with east and north
                if esist_nordest == 0 and esist_nord == 1 and esist_est == 1:
                    if (np.linalg.norm(Vertex[ind_central-1,:]-Vertex[ind_est-1,:]) <= TX_X84 and 
                        np.linalg.norm(Vertex[ind_central-1,:]-Vertex[ind_nord-1,:]) <= TX_X84 and 
                        np.linalg.norm(Vertex[ind_nord-1,:]-Vertex[ind_est-1,:]) <= TX_X84):
                        ind += 1
                        Triangle[ind-1,:] = [ind_central-1, ind_est-1, ind_nord-1]  # Convert to 0-based indexing
    
    Triangle = Triangle[:ind,:]
    return Triangle, Vertex

