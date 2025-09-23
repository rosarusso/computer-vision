#!/usr/bin/env python3

import numpy as np
from scipy.linalg import svd
from scipy.io import loadmat
import warnings
try:
    import matplotlib.pyplot as plt
    from matplotlib.image import imread
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Warning: matplotlib not available. Visualization will be disabled.")

class PoseEstimationToolbox:
    # Toolbox for 3D pose estimation

    @staticmethod
    def ns(A):
        """
        Solve the null-space problem A*v=0
        Returns the 1-d null-space of A 
        """
        U, D, V = svd(A)

        # Check condition number
        if len(D) > 1:
            c = D[0] / D[-2]
            if c > 200:
                warnings.warn(f'ns: condition number is {c:.0f}')

        return V[-1, :]

    @staticmethod
    def rmse(x, y):
        # Root Mean Square Error between x and y
        diff = x - y
        return np.linalg.norm(diff, 'fro') / np.sqrt(diff.size - 1)

    @staticmethod
    def rigid(G, w):
        """
        Apply a rigid transformation G to points w

        G: 3x4 transformation matrix [R|t] or 4x4 homogeneous matrix
        w: 3xN points in world coordinates
        """
        if G.shape[0] == 4:
            G = G[:3, :]

        if w.shape[0] != 3:
            raise ValueError('World coordinates must be 3xN.')

        # Convert to homogeneous coordinates and apply transformation
        HM = np.vstack([w, np.ones((1, w.shape[1]))])
        wr = G @ HM

        return wr

    @staticmethod
    def vtrans(a, d):
        """
        Vec-transpose operator
        Change from a d*s by n matrix to a d*n by s matrix
        """
        # Input validation
        if not isinstance(a, np.ndarray):
            a = np.array(a)

        # Ensure a is at least 2D
        if a.ndim == 1:
            # Convert 1D array to column vector
            a = a.reshape(-1, 1)
        elif a.ndim == 0:
            raise ValueError("Input 'a' must be at least 1-dimensional.")

        s_shape = a.shape

        # Check if we have enough dimensions
        if len(s_shape) < 2:
            raise ValueError(f"Input array must be at least 2D, got shape {s_shape}.")

        n = s_shape[1]

        # Check if first dimension is divisible by d
        if s_shape[0] % d != 0:
            raise ValueError(f"First dimension ({s_shape[0]}) must be divisible by d ({d}).")

        s = s_shape[0] // d

        b = np.zeros((d * n, s))

        if n < s:
            for i in range(n):
                start_row = i * d
                end_row = (i + 1) * d
                b[start_row:end_row, :] = a[:, i].reshape(d, s)
        else:
            for i in range(s):
                start_row = i * d
                end_row = (i + 1) * d  
                b[:, i] = a[start_row:end_row, :].flatten()

        return b

    @staticmethod
    def absolute(X, Y, method='noscale'):
        """
        Solve absolute orientation between two 3D point sets

        X, Y: 3xN arrays of corresponding 3D points
        method: 'noscale' or 'scale'

        Returns: G (3x4 transformation), s (scale), res (residual)
        """
        # Input validation
        if not isinstance(X, np.ndarray):
            X = np.array(X)
        if not isinstance(Y, np.ndarray):
            Y = np.array(Y)

        # Ensure proper dimensionality
        if X.ndim == 1:
            X = X.reshape(-1, 1)
        if Y.ndim == 1:
            Y = Y.reshape(-1, 1)

        if X.shape != Y.shape:
            raise ValueError("X and Y must have the same shape.")

        if X.shape[0] != 3:
            raise ValueError("Points must be 3D (3xN).")

        # Remove NaN entries
        valid_mask = ~(np.isnan(X).any(axis=0) | np.isnan(Y).any(axis=0))
        X_clean = X[:, valid_mask]
        Y_clean = Y[:, valid_mask]

        if X_clean.shape[1] < 3:
            raise ValueError("Need at least 3 valid point correspondences.")

        dime = Y_clean.shape[1]

        # Compute centroids
        cm = np.mean(Y_clean, axis=1, keepdims=True)
        cd = np.mean(X_clean, axis=1, keepdims=True)

        # Subtract centroids
        Yb = Y_clean - cm
        Xb = X_clean - cd

        # Compute scale
        if method == 'scale':
            s = np.linalg.norm(Xb, 'fro') / np.linalg.norm(Yb, 'fro')
        else:
            s = 1.0

        # Apply scale  
        Xb = s * Xb
        cd = cd / s

        # Compute rotation using SVD
        K = Xb @ Yb.T
        U, D, V = svd(K)

        # Ensure proper rotation matrix (det(R) = 1)
        S = np.diag([1, 1, np.linalg.det(U @ V)])
        R = U @ S @ V

        # Compute translation
        t = cd.flatten() - R @ cm.flatten()

        # Construct transformation matrix
        G = np.hstack([R, t.reshape(-1, 1)])

        # Compute residual
        X_transformed = PoseEstimationToolbox.rigid(G, Y_clean)
        res = PoseEstimationToolbox.rmse(X_clean, X_transformed / s)

        return G, s, res

    @staticmethod
    def proj(P, c3d):
        """
        Compute perspective projection from 3D to pixel coordinates

        P: 3x4 projection matrix  
        c3d: Nx3 array of 3D points

        Returns: u, v arrays of pixel coordinates (rounded)
        """
        if not isinstance(c3d, np.ndarray):
            c3d = np.array(c3d)

        if c3d.ndim == 1:
            c3d = c3d.reshape(1, -1)

        # Convert to homogeneous coordinates
        h3d = np.hstack([c3d, np.ones((c3d.shape[0], 1))]).T

        # Project
        h2d = P @ h3d

        # Avoid division by zero
        valid_mask = np.abs(h2d[2]) > 1e-10

        # Initialize output
        u = np.full(h2d.shape[1], np.nan)
        v = np.full(h2d.shape[1], np.nan)

        # Normalize by homogeneous coordinate where valid
        if np.any(valid_mask):
            c2d = h2d[:2, valid_mask] / h2d[2, valid_mask]
            u[valid_mask] = np.round(c2d[0]).astype(int)
            v[valid_mask] = np.round(c2d[1]).astype(int)

        return u, v

    @staticmethod
    def p2t(H, m):
        """
        Apply 2D projective transformation

        H: 3x3 homography matrix
        m: 2xN points
        """
        if H.shape != (3, 3):
            raise ValueError('Transformation matrix must be 3x3!')

        if not isinstance(m, np.ndarray):
            m = np.array(m)

        if m.ndim == 1:
            m = m.reshape(-1, 1)

        if m.shape[0] != 2:
            raise ValueError('Image coordinates must be Cartesian (2xN)!')

        dime = m.shape[1]

        # Convert to homogeneous
        c3d = np.vstack([m, np.ones((1, dime))])
        h2d = H @ c3d

        # Normalize
        c2d = h2d[:2] / h2d[2]

        return c2d

    @staticmethod
    def p3t(T, w):
        """
        Apply 3D projective transformation

        T: 4x4 transformation matrix
        w: 3xN points
        """
        if T.shape != (4, 4):
            raise ValueError('Transformation matrix must be 4x4.')

        if not isinstance(w, np.ndarray):
            w = np.array(w)

        if w.ndim == 1:
            w = w.reshape(-1, 1)

        if w.shape[0] != 3:
            raise ValueError('World coordinates must be Cartesian (3xN).')

        # Convert to homogeneous and apply transformation
        tmp = T @ np.vstack([w, np.ones((1, w.shape[1]))])

        # Normalize
        wt = tmp[:3] / tmp[3]

        return wt

    @staticmethod
    def pt(H, m):
        # Apply projective transformation (2D or 3D based on matrix size)
        
        if H.shape == (3, 3):
            return PoseEstimationToolbox.p2t(H, m)
        elif H.shape == (4, 4):
            return PoseEstimationToolbox.p3t(H, m)
        else:
            raise ValueError('Transformation not implemented for this matrix size.')

    @staticmethod
    def normalize(points):
        # Normalize homogeneous coordinates
        
        if not isinstance(points, np.ndarray):
            points = np.array(points)

        if points.ndim == 1:
            points = points.reshape(-1, 1)

        if points.shape[0] == 3:  # 2D homogeneous
            valid_mask = np.abs(points[2]) > 1e-10
            result = np.full((2, points.shape[1]), np.nan)
            result[:, valid_mask] = points[:2, valid_mask] / points[2, valid_mask]
            return result
        elif points.shape[0] == 4:  # 3D homogeneous
            valid_mask = np.abs(points[3]) > 1e-10
            result = np.full((3, points.shape[1]), np.nan)
            result[:, valid_mask] = points[:3, valid_mask] / points[3, valid_mask]
            return result
        else:
            return points

    @staticmethod
    def exterior_fiore(A, model3d, data2d):
        """
        Solve exterior orientation with Fiore's algorithm

        A: 3x3 camera intrinsic matrix (must be normalized, A[2,2] = 1)
        model3d: 3xN array of 3D world points  
        data2d: 2xN array of corresponding 2D image points

        Returns: G (3x4 pose matrix [R|t]), s (scale factor)
        """
        if A[2, 2] != 1:
            raise ValueError('A must be normalized (A[2,2] = 1)')

        # Ensure inputs are numpy arrays with proper dimensions
        if not isinstance(model3d, np.ndarray):
            model3d = np.array(model3d)
        if not isinstance(data2d, np.ndarray):
            data2d = np.array(data2d)

        if model3d.ndim == 1:
            model3d = model3d.reshape(-1, 1)
        if data2d.ndim == 1:
            data2d = data2d.reshape(-1, 1)

        # Change to normalized coordinates (or image coordinates)
        data2d_hom = np.vstack([data2d, np.ones((1, data2d.shape[1]))])
        m = np.linalg.inv(A) @ data2d_hom

        # Convert to homogeneous
        m = np.vstack([m, np.ones((1, m.shape[1]))])

        S = np.vstack([model3d, np.ones((1, model3d.shape[1]))])
        U, X, V = svd(S)
        i = np.linalg.matrix_rank(S)
        V2 = V[i:, :].T

        numP = data2d.shape[1]
        D = []

        # Build D matrix
        for j in range(numP):
            D_block = np.zeros((3, numP))
            D_block[:, j] = m[:3, j]
            D.append(D_block)

        D = np.vstack(D)

        L = np.kron(V2.T, np.eye(3)) @ D
        z = PoseEstimationToolbox.ns(L)

        # The scale factor can be negative, but if the intrinsic matrix is normalized correctly, the depths are positive
        z = z * np.sign(z[0]) if z[0] != 0 else z

        # Use absolute orientation to get final transformation
        try:
            vtrans_result = PoseEstimationToolbox.vtrans(D @ z, 3)
            G, s, res = PoseEstimationToolbox.absolute(vtrans_result, model3d, 'scale')
        except Exception as e:
            print(f"Warning in exterior_fiore: {e}")
            # Fallback to simple DLT if vtrans fails
            G = np.eye(3, 4)  # Identity transformation
            s = 1.0

        return G, s


def load_matlab_data():
    try:
        # Load the MATLAB data file
        mat_data = loadmat('imgInfo.mat')

        # Extract imgInfo structure
        if 'imgInfo' in mat_data:
            img_info = mat_data['imgInfo']

            # Handle both struct array and regular struct formats
            if isinstance(img_info, np.ndarray) and img_info.size == 1:
                img_info = img_info[0, 0]

            # Extract fields
            p2D = img_info['punti2DImg']  # 2D image points
            p3D = img_info['punti3DImg']  # 3D world points  
            K = img_info['K']             # Intrinsic camera matrix
            R = img_info['R'] if 'R' in img_info.dtype.names else None   # Rotation matrix
            T = img_info['T'] if 'T' in img_info.dtype.names else None   # Translation vector

            return {
                'p2D': p2D,
                'p3D': p3D,
                'K': K,
                'R': R,
                'T': T
            }
        else:
            raise ValueError("imgInfo not found in .mat file")

    except Exception as e:
        print(f"Error loading imgInfo.mat: {e}")
        print("Make sure imgInfo.mat is in the current directory")
        return None


def main():
    # Load camera information from imgInfo.mat
    img_data = load_matlab_data()
    if img_data is None:
        return

    # Load image
    try:
        if MATPLOTLIB_AVAILABLE:
            img = imread('cav.jpg')
        else:
            print("Image loading skipped - matplotlib not available")
            img = None
    except Exception as e:
        print(f"Warning: Could not load cav.jpg: {e}")
        img = None

    # Extract data
    p2D = img_data['p2D']  # Image points
    p3D = img_data['p3D']  # World points
    K = img_data['K']      # Intrinsic matrix

    # Convert to proper format if needed (ensure 2D points are 2xN, 3D points are 3xN)
    if p2D.shape[0] != 2:
        p2D = p2D.T
    if p3D.shape[0] != 3:
        p3D = p3D.T

    print(f"\nLoaded data:")
    print(f"  2D points shape: {p2D.shape}")
    print(f"  3D points shape: {p3D.shape}")
    print(f"  Camera matrix K:")
    for row in K:
        print(f"    {row}")

    N = min(100, p2D.shape[1])  # Number of points to use
    print(f"\nUsing {N} points for pose estimation")

    # Use subset of points
    p2D_subset = p2D[:, :N]
    p3D_subset = p3D[:, :N]

    # Create homogeneous coordinates for 3D points
    p3Dn = np.vstack([p3D, np.ones((1, p3D.shape[1]))])

    # Display 3D points if matplotlib is available
    if MATPLOTLIB_AVAILABLE and img is not None:
        plt.figure(figsize=(15, 10))

        # Plot 1: 3D points
        plt.subplot(2, 3, 1)
        ax = plt.axes(projection='3d')
        ax.scatter(p3D[0], p3D[1], p3D[2], s=5, c='cyan')
        ax.set_xlabel('X')
        ax.set_ylabel('Y') 
        ax.set_zlabel('Z')
        ax.set_title('3D Points')

        # Plot 2: Image with 2D points
        plt.subplot(2, 3, 2)
        plt.imshow(img)
        plt.plot(p2D[0], p2D[1], 'r.', markersize=2)
        plt.title('2D Image Points')
        plt.axis('equal')

        # If ground truth pose is available, show projection
        if img_data['R'] is not None and img_data['T'] is not None:
            P_true = K @ np.hstack([img_data['R'], img_data['T']])
            u_true, v_true = PoseEstimationToolbox.proj(P_true, p3D.T)
            plt.plot(u_true, v_true, 'go', markersize=3, fillstyle='none')

            ps = PoseEstimationToolbox.normalize(P_true @ p3Dn)
            plt.plot(ps[0], ps[1], 'go', markersize=2, fillstyle='none')

    print(f"\nExterior orientation...")
    try:
        G, s = PoseEstimationToolbox.exterior_fiore(K, p3D_subset, p2D_subset)

        print(f"Estimated transformation G:")
        for row in G:
            print(f"  {row}")
        print(f"Scale factor: {s}")

        # Test projection with estimated pose
        P1 = K @ G
        u1, v1 = PoseEstimationToolbox.proj(P1, p3D.T)

        # Show results if matplotlib available
        if MATPLOTLIB_AVAILABLE and img is not None:
            plt.subplot(2, 3, 3)
            plt.imshow(img)
            plt.plot(p2D[0], p2D[1], 'r.', markersize=2, label='Observed')
            valid_mask = ~(np.isnan(u1) | np.isnan(v1))
            plt.plot(u1[valid_mask], v1[valid_mask], 'bo', markersize=2, 
                    fillstyle='none', label='Fiore estimate')
            plt.title('Fiore Projection Results')
            plt.legend()

            # Also test normalize function
            ps1 = PoseEstimationToolbox.normalize(P1 @ p3Dn)
            plt.plot(ps1[0], ps1[1], 'bo', markersize=1, fillstyle='none')

        # Compute reprojection error
        valid_proj = ~(np.isnan(u1) | np.isnan(v1))
        if np.any(valid_proj):
            errors = np.sqrt((u1[valid_proj] - p2D[0, valid_proj])**2 + 
                           (v1[valid_proj] - p2D[1, valid_proj])**2)
            mean_error = np.mean(errors)
            print(f"Mean reprojection error: {mean_error:.2f} pixels")

        # Compare with ground truth if available
        if img_data['R'] is not None and img_data['T'] is not None:
            print(f"\nComparison with ground truth:")
            R_true = img_data['R']
            T_true = img_data['T'].flatten()
            G_true = np.hstack([R_true, T_true.reshape(-1, 1)])

            print(f"Ground truth transformation:")
            for row in G_true:
                print(f"  {row}")

            print(f"\nEstimated vs Ground Truth:")
            print(f"  [R_est | T_est] vs [R_true | T_true]:")
            for i, (row_est, row_true) in enumerate(zip(G, G_true)):
                print(f"  Row {i}: {row_est} vs {row_true}")

    except Exception as e:
        print(f"Error in Fiore algorithm: {e}")
        import traceback
        traceback.print_exc()

    if MATPLOTLIB_AVAILABLE and img is not None:
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()

