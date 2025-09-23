import numpy as np
import scipy.io as sio
from scipy.ndimage import map_coordinates
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import matplotlib
matplotlib.use('TkAgg')

def vec(A):
    # Convert matrix to column vector by stacking columns.
    return A.flatten(order='F').reshape(-1, 1)

def ivec(v, r):
    # Reshape column vector into matrix with r rows.
    v = np.array(v).flatten()
    if len(v) % r != 0:
        raise ValueError("Number of rows must divide length of vector.")
    return v.reshape(r, -1, order='F')

def p2t(H, m):
    # Apply 2D projective transformation (homography) to 2xN points.
    H = np.array(H)
    m = np.array(m)
    if H.shape != (3, 3):
        raise ValueError("Transformation matrix must be 3x3.")
    if m.shape[0] != 2:
        raise ValueError("Image coordinates must be Cartesian (2xN).")
    
    N = m.shape[1]
    # Convert to homogeneous
    c3d = np.vstack([m, np.ones((1, N))])  # 3xN
    # Apply homography
    h2d = H @ c3d  # 3xN
    # Convert back to Cartesian
    c2d = h2d[:2, :] / h2d[2, :]  # Divide by last coordinate
    return c2d

def ns(A):
    # Compute null space of matrix A using SVD.
    _, _, Vt = np.linalg.svd(A)
    return Vt[-1, :].reshape(-1, 1)

def KRt_fact(P):
    # Factorize P = K @ R @ [I | t]
    Q = P[:, :3]
    q = P[:, 3]
    
    # QR decomposition of inv(Q)
    O, T = np.linalg.qr(np.linalg.inv(Q))
    
    # Internal calibration matrix
    K = np.linalg.inv(T)
    
    # Sign correction for rotation
    detO = np.linalg.det(O)
    R = detO * O.T
    
    # Translation
    t = T @ (detO * q)
    
    # Normalize K so that K[2,2] = 1
    K = K / K[2, 2]
    
    return K, R, t

def fundamental(pl, pr):
    # Compute fundamental matrix and epipoles from two camera matrices.
    # Extract camera centers
    Ql = pl[:, :3]
    ql = pl[:, 3]
    Qr = pr[:, :3]
    qr = pr[:, 3]
    
    # Camera centers in world coordinates
    cl = -np.linalg.inv(Ql) @ ql
    cr = -np.linalg.inv(Qr) @ qr
    
    # Epipoles (projected camera centers in the other image)
    el = pl @ np.hstack([cr, 1])  # Left epipole: projection of right camera center in left image
    er = pr @ np.hstack([cl, 1])  # Right epipole: projection of left camera center in right image
    
    # Normalize epipoles
    el = el[:2] / el[2]
    er = er[:2] / er[2]
    
    # Skew-symmetric matrix of left epipole
    Sel = np.array([[ 0,     -1,     er[1]],
                    [ 1,      0,    -er[0]],
                    [-er[1], er[0],    0   ]])
    
    # Fundamental matrix
    F = Sel @ pr @ np.linalg.pinv(pl)
    
    return F, el, er


def imwarp(I, H, meth='linear', sz='same'):
    """
    Warp image I using homography H.
    
    Parameters:
        I: input image (H x W or H x W x C)
        H: 3x3 homography matrix
        meth: interpolation method ('linear' only supported)
        sz: 'same', 'valid', or [minx, miny, maxx, maxy]
    
    Returns:
        I2: warped image
        bb: bounding box [minx; miny; maxx; maxy]
        alpha: mask of valid pixels (if requested)
    """
    H = np.array(H)
    if H.shape != (3, 3):
        raise ValueError("Invalid transformation matrix: must be 3x3")
    
    I = np.atleast_3d(I).astype(float)
    h, w = I.shape[:2]
    
    # Default parameters
    if meth not in ['linear']:
        raise ValueError("Only 'linear' interpolation is supported.")
    
    # Determine bounding box
    if sz == 'same':
        minx, maxx = 0, w - 1
        miny, maxy = 0, h - 1
    elif sz == 'valid':
        corners = np.array([[0, 0, w, w],
                            [0, h, 0, h]])  # 2x4
        corners_x = p2t(H, corners)
        minx = np.floor(corners_x[0, :].min()).astype(int)
        maxx = np.ceil(corners_x[0, :].max()).astype(int)
        miny = np.floor(corners_x[1, :].min()).astype(int)
        maxy = np.ceil(corners_x[1, :].max()).astype(int)
    elif isinstance(sz, (list, tuple, np.ndarray)) and len(sz) == 4:
        minx, miny, maxx, maxy = map(int, sz)
    else:
        raise ValueError("Invalid size option.")
    
    bb = np.array([minx, miny, maxx, maxy])
    
    # Create mesh grid in output space
    x = np.arange(minx, maxx)
    y = np.arange(miny, maxy)
    X, Y = np.meshgrid(x, y)
    
    # Flatten coordinates
    flat_X = X.flatten()
    flat_Y = Y.flatten()
    
    # Stack as 2xN
    out_coords = np.vstack([flat_X, flat_Y])  # 2xN
    
    # Apply inverse homography to get source coordinates
    H_inv = np.linalg.inv(H)
    src_coords_homo = p2t(H_inv, out_coords)  # 2xN
    
    # Prepare for map_coordinates
    src_x = src_coords_homo[0, :]
    src_y = src_coords_homo[1, :]
    
    # Interpolate
    I2 = np.zeros((Y.shape[0], X.shape[1], I.shape[2]), dtype=float)
    alpha = np.zeros_like(I2[:, :, 0], dtype=bool)
    
    for c in range(I.shape[2]):
        # Interpolate using map_coordinates
        values = map_coordinates(I[:, :, c], [src_y, src_x], order=1, mode='constant', cval=np.nan, prefilter=False)
        values = values.reshape(Y.shape)
        I2[:, :, c] = values
        alpha |= (~np.isnan(values))
    
    # Convert to original dtype
    I2 = np.where(np.isnan(I2), 0, I2)
    I2 = I2.squeeze()
    
    return I2, bb, alpha

def main():
    # Load calibration data
    stereo = sio.loadmat('Calib_Results_stereo.mat')
    
    # Load images
    imgleft = plt.imread('left.jpg')
    imgright = plt.imread('right.jpg')
    
    # Handle grayscale images
    if imgleft.ndim == 2:
        imgleft = np.atleast_3d(imgleft)
    if imgright.ndim == 2:
        imgright = np.atleast_3d(imgright)
    
    # Ensure images are float for processing
    imgleft = imgleft.astype(float)
    imgright = imgright.astype(float)
    
    # Extract calibration matrices
    KK_left = stereo['KK_left']
    KK_right = stereo['KK_right']
    R_stereo = stereo['R']
    T_stereo = stereo['T'].flatten()
    
    # Construct projection matrices
    P_left = KK_left @ np.hstack([np.eye(3), np.zeros((3, 1))])  # K @ [I | 0]
    P_right = KK_right @ np.hstack([R_stereo, T_stereo.reshape(-1, 1)])  # K @ [R | t]
    
    Q_left = P_left[:, :3]
    Q_right = P_right[:, :3]
    q_left = P_left[:, 3]
    q_right = P_right[:, 3]
    
    # Factorize P matrices
    Kl, R_left, t_left = KRt_fact(P_left)
    Kr, R_right, t_right = KRt_fact(P_right)
    
    # Camera centers
    c_left = -np.linalg.inv(Q_left) @ q_left
    c_right = -np.linalg.inv(Q_right) @ q_right
    
    # Baseline vector
    v0 = c_right - c_left
    
    # New rotation matrix
    k = R_left[2, :]  # third row
    v1 = np.cross(k, v0)
    v2 = np.cross(v0, v1)
    
    # Normalize and build R
    R = np.array([
        v0 / np.linalg.norm(v0),
        v1 / np.linalg.norm(v1),
        v2 / np.linalg.norm(v2)
    ])
    
    # New projection matrices (rectified)
    Pnleft = Kl @ np.hstack([R, -R @ c_left.reshape(3, 1)])
    Pnright = Kr @ np.hstack([R, -R @ c_right.reshape(3, 1)])
    
    # Homographies (warping matrices)
    T_left = Pnleft[:, :3] @ np.linalg.inv(Q_left)
    T_right = Pnright[:, :3] @ np.linalg.inv(Q_right)
    
    # Warp images
    imgleft_rect, bb_left, _ = imwarp(imgleft, T_left, meth='linear', sz='same')
    imgright_rect, bb_right, _ = imwarp(imgright, T_right, meth='linear', sz='same')
    imgright_rect = imgright_rect[::-1, ::-1]
    
    # Compute fundamental matrix
    F, eleft, eright = fundamental(Pnleft, Pnright)
    
    # Display images
    fig1 = plt.figure(1)
    plt.show(block=False)
    plt.imshow(np.clip(imgleft_rect, 0, 255).astype(np.uint8))
    plt.title('Click 3 points on left image')
    plt.axis('image')
    
    fig2 = plt.figure(2)
    plt.xlim(0, imgright_rect.shape[1])
    plt.ylim(imgright_rect.shape[0], 0)
    plt.show(block=False)
    plt.imshow(np.clip(imgright_rect, 0, 255).astype(np.uint8))
    plt.title('Corresponding Epipolar Lines')
    plt.axis('image')
    
    colors = ['yellow', 'green', [1, 0.6, 0.1]]
    
    # Click 3 points
    for i in range(3):
        print(f"Click point {i+1}")
        plt.figure(1)
        pt = plt.ginput(1, timeout=30)
        if not pt:
            print("Timeout or no input.")
            break
        x_left, y_left = pt[0]
        
        # Plot point on left image
        plt.plot(x_left, y_left, 'rx', markersize=10)
        plt.draw()
        
        # Compute epipolar line in right image: l' = F @ p
        p_left = np.array([x_left, y_left, 1])
        line_right = F @ p_left  # [a, b, c]

        # Extract line parameters: ax + by + c = 0
        a, b, c = line_right

        # Get image dimensions
        h, w = imgright_rect.shape[:2]

        # Handle degenerate cases
        if abs(b) < 1e-6 and abs(a) < 1e-6:
            continue  # Invalid line

        # Calculate intersection points with image borders
        x_vals = []
        y_vals = []

        # Left border (x = 0)
        if abs(b) > 1e-6:
            y = -c / b
            if 0 <= y < h:
                x_vals.append(0)
                y_vals.append(y)

        # Right border (x = w - 1)
        if abs(b) > 1e-6:
            y = (-a * (w - 1) - c) / b
            if 0 <= y < h:
                x_vals.append(w - 1)
                y_vals.append(y)

        # Top border (y = 0)
        if abs(a) > 1e-6:
            x = -c / a
            if 0 <= x < w:
                x_vals.append(x)
                y_vals.append(0)

        # Bottom border (y = h - 1)
        if abs(a) > 1e-6:
            x = (-b * (h - 1) - c) / a
            if 0 <= x < w:
                x_vals.append(x)
                y_vals.append((h - 1))
        # Only plot if we have valid points
        if len(x_vals) >= 2:
            plt.figure(2)
            y_vals = [imgright_rect.shape[0] - k for k in y_vals]
            plt.plot(x_vals[:2], y_vals[:2], color=colors[i], linewidth=2)
            plt.draw()


    plt.show()

if __name__ == "__main__":
    main()

