import numpy as np
import matplotlib.pyplot as plt
import cv2
from scipy.linalg import svd
from scipy.io import loadmat

def ns(A):
    # Compute nullspace of matrix A using SVD
	
    _, _, V = svd(A)
    return V[:, -1]

def intersect_base(PPM, m):
    # Triangulate 3D points from 2D image correspondences
    numP = len(m[0])
    numV = len(m)
    
    M = []
    for i in range(numP):
        A = []
        for view in range(numV):
            # Normalize projection matrix
            PPM[view] = PPM[view] / np.linalg.norm(PPM[view][2, :3])
            
            # Build constraint equations
            row1 = PPM[view][0, :] - m[view][i][0] * PPM[view][2, :]
            row2 = PPM[view][1, :] - m[view][i][1] * PPM[view][2, :]
            A.append(row1)
            A.append(row2)
        
        # Stack constraints into matrix
        A = np.array(A)
        
        # Compute nullspace solution
        x = ns(A)
        
        # Convert to inhomogeneous coordinates
        M.append(x[:3] / x[3])
    
    return np.array(M).T

# Load projection matrices from .mat files
calib1 = loadmat('Calib_direct_1.mat')
P1 = calib1['P']  # Assuming the variable is named 'P' in the .mat file

calib2 = loadmat('Calib_direct_2.mat')
P2 = calib2['P']  # Assuming the variable is named 'P' in the .mat file

# Projection matrices
PPM = [P1, P2]

# Load images
I1 = cv2.imread('Image1.jpg')
I2 = cv2.imread('Image2.jpg')

# Convert BGR to RGB for matplotlib
I1 = cv2.cvtColor(I1, cv2.COLOR_BGR2RGB)
I2 = cv2.cvtColor(I2, cv2.COLOR_BGR2RGB)

# Pick points in first image
plt.figure()
plt.imshow(I1)
plt.title('Pick 2 points')
p1 = plt.ginput(1)[0]
plt.plot(p1[0], p1[1], 'g*')
p2 = plt.ginput(1)[0]
plt.plot(p2[0], p2[1], 'g*')
mleft = [p1, p2]
plt.show()

# Pick corresponding points in second image
plt.figure()
plt.imshow(I2)
plt.title('Pick the same points of the previous image')
p1 = plt.ginput(1)[0]
plt.plot(p1[0], p1[1], 'g*')
p2 = plt.ginput(1)[0]
plt.plot(p2[0], p2[1], 'g*')
mright = [p1, p2]
plt.show()

# Prepare correspondences
m = [mleft, mright]

# Triangulate points
M = intersect_base(PPM, m)

# Calculate Euclidean distance
d = np.linalg.norm(M[:, 0] - M[:, 1])

# Display results
plt.figure()
plt.imshow(I1)
# Plot selected points
x_coords = [point[0] for point in mleft]
y_coords = [point[1] for point in mleft]
plt.plot(x_coords, y_coords, 'g*')
# Add distance text
mid_x = (mleft[0][0] + mleft[1][0]) / 2 + 10
mid_y = (mleft[0][1] + mleft[1][1]) / 2 + 10
plt.text(mid_x, mid_y, f'{d:.2f}', fontsize=16, color='green')
plt.show()

print(f"3D Points:\n{M}")
print(f"Distance between points: {d:.2f}")

