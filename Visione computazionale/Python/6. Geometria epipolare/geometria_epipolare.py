import numpy as np
import matplotlib.pyplot as plt
import cv2
from scipy.io import loadmat
from scipy.linalg import qr, svd

def fundamental(pl, pr):
    # Compute fundamental matrix from projection matrices

    # Optical centers
    cl = -np.linalg.inv(pl[:, :3]) @ pl[:, 3]
    cr = -np.linalg.inv(pr[:, :3]) @ pr[:, 3]
    
    # Epipolar lines
    el = pl @ np.append(cr, 1)
    er = pr @ np.append(cl, 1)
    
    el = el / np.linalg.norm(el)
    er = er / np.linalg.norm(er)
    
    # Fundamental Matrix
    cross_er = np.array([
        [0, -er[2], er[1]],
        [er[2], 0, -er[0]],
        [-er[1], er[0], 0]
    ])
    
    F = cross_er @ pr[:, :3] @ np.linalg.inv(pl[:, :3])
    F = F / np.linalg.norm(F)
    
    return F, el, er

def ns(A):
    #Compute nullspace of matrix A
    
    _, _, V = svd(A)
    return V[:, -1]

def intersect_base(PPM, m):
    # Triangulate points using multiple views
    
    numP = len(m[0])
    numV = len(m)
    
    M = []
    for i in range(numP):
        A = []
        for view in range(numV):
            P = PPM[view]
            # Normalize P
            P = P / np.linalg.norm(P[2, :3])
            A.append(P[0, :] - m[view][i][0] * P[2, :])
            A.append(P[1, :] - m[view][i][1] * P[2, :])
        
        A = np.array(A)
        x = ns(A)
        M.append(x[:3] / x[3])
    
    return np.array(M).T

def KRt_fact(P):
    # Factorize projection matrix P into K, R, and t
    
    Q = P[:, :3]
    q = P[:, 3]
    
    O, T = qr(np.linalg.inv(Q))
    
    # Internal parameters
    K = np.linalg.inv(T)
    detO = np.linalg.det(O)
    
    # Sign correction
    R = detO * O.T
    
    # Translation
    t = T @ (detO * q)
    
    # Normalization
    K = K / K[2, 2]
    
    return K, R, t

def vec(A):
    # Convert matrix to column vector
   
    return A.reshape(-1, 1)

def ivec(v, r):
    # Convert column vector to matrix with r rows
    
    if v.ndim > 1 and v.shape[1] != 1:
        raise ValueError("Input vector must be a column vector!")
    
    if len(v) % r != 0:
        raise ValueError("Number of rows is not compatible with vector length!")
    
    return v.reshape(r, -1)

def p2t(H, m):
    # Apply 2D projective transformation
    
    if H.shape != (3, 3):
        raise ValueError("Transformation matrix must be 3x3!")
    
    if m.shape[0] != 2:
        raise ValueError("Image coordinates must be Cartesian (2 rows)!")
    
    dime = m.shape[1]
    c3d = np.vstack([m, np.ones((1, dime))])
    h2d = H @ c3d
    c2d = h2d[:2, :] / h2d[2, :]
    
    return c2d

def load_ppm(filename):
    # Load PPM file
    
    data = loadmat(filename)
    return data['P']

def draw_epipolar_lines():
    # Load projection matrices
    try:
        pl = load_ppm('Data/angel1.ppm')
        pr = load_ppm('Data/angel2.ppm')
    except:
        print("PPM files not found. Using identity matrices for demonstration.")
        pl = np.eye(3, 4)
        pr = np.eye(3, 4)
        pr[0, 3] = 0.5  # Add some translation for variety
    
    # Compute fundamental matrix
    F, ep_left, ep_right = fundamental(pl, pr)
    
    # Load images
    try:
        imgleft = cv2.imread('Data/angel1.JPG')
        imgright = cv2.imread('Data/angel2.JPG')
        if imgleft is None or imgright is None:
            raise FileNotFoundError
    except:
        print("Image files not found.")
    
    imgleft = cv2.cvtColor(imgleft, cv2.COLOR_BGR2RGB)
    imgright = cv2.cvtColor(imgright, cv2.COLOR_BGR2RGB)
    
    # Create figures
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    ax1.imshow(imgleft)
    ax1.set_title('Click 3 points')
    ax1.axis('image')
    
    ax2.imshow(imgright)
    ax2.set_title('Corresponding Epipolar Lines')
    ax2.axis('image')
    
    colors = ['yellow', 'green', 'red']
    
    # Point picking and epipolar line calculation
    points = []
    print("Click 3 points on the image")
    
    
    print("Click 3 points (the same as earlier).")
    points = plt.ginput(3)
    points = [p for p in points]  # Convert to list of tuples

    # Now plot all points and epipolar lines after collection
    for i, point in enumerate(points):
        # Plot point on left image
        ax1.plot(point[0], point[1], 'rx')

        # Calculate epipolar line on right image
        p_left = np.array([point[0], point[1], 1])
        p_right = F @ p_left

        # Generate epipolar line
        n = imgright.shape[1]  # Use imgright width
        epipolar_x = np.array([0, n - 1])
        epipolar_y = (-p_right[2] - p_right[0] * epipolar_x) / (p_right[1] + 1e-6)  # Avoid div by zero

        # Plot epipolar line on right image
        ax2.plot(epipolar_x, epipolar_y, color=colors[i], linewidth=1)

    plt.draw()   # Force redraw
    plt.show()   # Show final result
    return points

    return points

def triangulate_points():
    
    # Load projection matrices
    try:
        data1 = loadmat('Calib_direct_1.mat')
        P1 = data1['P']
        data2 = loadmat('Calib_direct_2.mat')
        P2 = data2['P']
    except:
        print("Calibration files not found.")
    
    # Projection Matrix collection
    PPM = [P1, P2]
    
    # Load images
    try:
        I1 = cv2.imread('Image1.jpg')
        I2 = cv2.imread('Image2.jpg')
        if I1 is None or I2 is None:
            raise FileNotFoundError
    except:
        print("Image files not found.")
        
    I1 = cv2.cvtColor(I1, cv2.COLOR_BGR2RGB)
    I2 = cv2.cvtColor(I2, cv2.COLOR_BGR2RGB)
    
    # Pick points on first image
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    ax1.imshow(I1)
    ax1.set_title('Pick 2 points on Image 1')
    
    mleft = []
    for i in range(2):
        point = plt.ginput(1)[0]
        mleft.append(point)
        ax1.plot(point[0], point[1], 'g*')
    
    plt.show()
    
    # Pick corresponding points on second image
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    ax2.imshow(I2)
    ax2.set_title('Pick corresponding 2 points on Image 2')
    
    mright = []
    for i in range(2):
        point = plt.ginput(1)[0]
        mright.append(point)
        ax2.plot(point[0], point[1], 'g*')
    
    plt.show()
    
    # Group measurements
    m = [mleft, mright]
    
    # Model estimation
    M = intersect_base(PPM, m)
    
    # Calculate Euclidean distance between points
    d = np.linalg.norm(M[:, 0] - M[:, 1])
    
    # Display results
    fig3, ax3 = plt.subplots(figsize=(8, 6))
    ax3.imshow(I1)
    ax3.set_title(f'Distance: {d:.2f}')
    
    # Draw lines between points
    x_coords = [point[0] for point in mleft]
    y_coords = [point[1] for point in mleft]
    ax3.plot(x_coords, y_coords, 'g-', linewidth=2)
    ax3.plot(x_coords, y_coords, 'go')
    
    # Add distance text
    mid_x = (mleft[0][0] + mleft[1][0]) / 2 + 10
    mid_y = (mleft[0][1] + mleft[1][1]) / 2 + 10
    ax3.text(mid_x, mid_y, f'{d:.2f}', fontsize=16, color='green')
    
    plt.show()
    
    return d, M

def main():
    points = draw_epipolar_lines()

if __name__ == "__main__":
    main()

