import numpy as np
import cv2
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

def vec(A):
    # Reorganize matrix elements into a column vector.
    return A.reshape(-1, 1, order='F')

def ivec(v, r):
    # Create a matrix grouping elements of a column vector.
    if len(v.shape) < 2 or v.shape[1] != 1:
        raise ValueError("Input vector must be a column vector.")
    
    length = v.shape[0]
    if length % r != 0:
        raise ValueError("Number of rows not suitable.")
    
    return v.reshape(r, length // r, order='F')

def p2t(H, m):
    # Apply a 2D projective transformation.
    if H.shape != (3, 3):
        raise ValueError("Transformation matrix must be 3x3")
    
    if m.shape[0] != 2:
        raise ValueError("Image coordinates must be Cartesian")
    
    dime = m.shape[1]
    c3d = np.vstack([m, np.ones((1, dime))])
    h2d = H @ c3d
    
    # Handle potential division by zero
    h2d[2, :] = np.where(np.abs(h2d[2, :]) < 1e-10, 1e-10, h2d[2, :])
    c2d = h2d[:2, :] / h2d[2, :]
    
    return c2d

def compute_homography_cv2(p1, p2):
    dst_pts = p1.T.astype(np.float32)  # Points from image 1
	src_pts = p2.T.astype(np.float32)  # Points from image 2
    
    # Compute homography from img2 to img1
    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0)
    
    if H is None:
        raise ValueError("Could not compute homography")
    
    return H

def compute_homography_dlt(p1, p2):
    n = p1.shape[1]
    if n < 4:
        raise ValueError("Need at least 4 point correspondences")
    
    # Normalize points for better numerical stability
    def normalize_points(pts):
        centroid = np.mean(pts, axis=1, keepdims=True)
        centered = pts - centroid
        scale = np.sqrt(2) / np.mean(np.sqrt(np.sum(centered**2, axis=0)))
        T = np.array([[scale, 0, -scale*centroid[0,0]], 
                      [0, scale, -scale*centroid[1,0]], 
                      [0, 0, 1]])
        return T, scale, centroid
    
    T1, scale1, centroid1 = normalize_points(p1)
    T2, scale2, centroid2 = normalize_points(p2)
    
    # Apply normalization
    p1_norm = T1 @ np.vstack([p1, np.ones((1, n))])
    p2_norm = T2 @ np.vstack([p2, np.ones((1, n))])
    p1_norm = p1_norm[:2, :]
    p2_norm = p2_norm[:2, :]
    
    A = []
    for i in range(n):
        x2, y2 = p2_norm[:, i]  # Source points (img2)
        x1, y1 = p1_norm[:, i]  # Destination points (img1)
        
        # Build constraint matrix for H such that H * p2 = p1
        A.append([0, 0, 0, -x2, -y2, -1, y1*x2, y1*y2, y1])
        A.append([x2, y2, 1, 0, 0, 0, -x1*x2, -x1*y2, -x1])
    
    A = np.array(A)
    
    # Solve using SVD
    _, _, Vt = np.linalg.svd(A)
    h = Vt[-1, :]  # Last row
    H_norm = h.reshape(3, 3)
    
    # Denormalize
    H = np.linalg.inv(T1) @ H_norm @ T2
    
    # Normalize so that H[2,2] = 1
    H = H / H[2, 2]
    
    return H

def imwarp_improved(I, H, method='linear'):
    h, w = I.shape
    
    # Find the bounding box of the warped image
    corners = np.array([[0, w-1, w-1, 0],
                       [0, 0, h-1, h-1]])
    
    warped_corners = p2t(H, corners)
    
    x_min = int(np.floor(np.min(warped_corners[0, :])))
    x_max = int(np.ceil(np.max(warped_corners[0, :])))
    y_min = int(np.floor(np.min(warped_corners[1, :])))
    y_max = int(np.ceil(np.max(warped_corners[1, :])))
    
    # Create output image
    out_h = y_max - y_min + 1
    out_w = x_max - x_min + 1
    
    # Create coordinate meshgrid for output image
    x_out, y_out = np.meshgrid(np.arange(x_min, x_max + 1),
                               np.arange(y_min, y_max + 1))
    
    # Map output coordinates back to input image
    coords_out = np.vstack([x_out.ravel(), y_out.ravel()])
    coords_in = p2t(np.linalg.inv(H), coords_out)
    
    x_in = coords_in[0, :].reshape(out_h, out_w)
    y_in = coords_in[1, :].reshape(out_h, out_w)
    
    # Create mask for valid coordinates
    valid = ((x_in >= 0) & (x_in < w-1) & (y_in >= 0) & (y_in < h-1))
    
    map_x = x_in.astype(np.float32)
    map_y = y_in.astype(np.float32)
    
    if method == 'linear':
        interpolation = cv2.INTER_LINEAR
    elif method == 'cubic':
        interpolation = cv2.INTER_CUBIC
    else:
        interpolation = cv2.INTER_LINEAR
    
    warped = cv2.remap(I.astype(np.float32), map_x, map_y, 
                       interpolation, borderMode=cv2.BORDER_CONSTANT, 
                       borderValue=0)
    
    # Apply validity mask
    warped = np.where(valid, warped, 0)
    
    bbox = [x_min, y_min, x_max, y_max]
    
    return warped, bbox, valid

def main():
    # Load images
    img1 = cv2.imread('city1.jpg', cv2.IMREAD_GRAYSCALE)
    img2 = cv2.imread('city2.jpg', cv2.IMREAD_GRAYSCALE)
    
    if img1 is None or img2 is None:
        raise FileNotFoundError("Could not load city1.jpg or city2.jpg")
    
    # Convert to float for better processing
    img1 = img1.astype(np.float64)
    img2 = img2.astype(np.float64)
    
    # Display images
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    ax1.imshow(img1, cmap='gray')
    ax1.set_title('Image 1 - Click 4 points')
    ax2.imshow(img2, cmap='gray')
    ax2.set_title('Image 2 - Click 4 corresponding points')
    plt.tight_layout()
    
    # Point picking (4 points minimum for homography)
    print("Please select 4 corresponding points in both images")
    print("Click on image 1 first, then corresponding points on image 2")
    points = plt.ginput(8, timeout=0)  # 8 total points (4 per image)
    plt.close()
    
    if len(points) != 8:
        raise ValueError("Need exactly 8 points (4 per image)")
    
    # Separate points for each image
    p1 = np.array([[p[0], p[1]] for p in points[:4]]).T  # Points in image 1
    p2 = np.array([[p[0], p[1]] for p in points[4:]]).T  # Corresponding points in image 2
    
    # Show selected points
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    ax1.imshow(img1, cmap='gray')
    ax1.scatter(p1[0, :], p1[1, :], c='red', marker='*', s=100)
    for i in range(p1.shape[1]):
        ax1.annotate(f'{i+1}', (p1[0, i], p1[1, i]), xytext=(5, 5), 
                    textcoords='offset points', color='red', fontweight='bold')
    ax1.set_title('Selected points on Image 1')
    
    ax2.imshow(img2, cmap='gray')
    ax2.scatter(p2[0, :], p2[1, :], c='green', marker='*', s=100)
    for i in range(p2.shape[1]):
        ax2.annotate(f'{i+1}', (p2[0, i], p2[1, i]), xytext=(5, 5), 
                    textcoords='offset points', color='green', fontweight='bold')
    ax2.set_title('Selected points on Image 2')
    plt.tight_layout()
    plt.show()
    
    try:
        # OpenCV method
        H = compute_homography_cv2(p1, p2)
        print("Using OpenCV homography computation")
    except:
        # Fallback to DLT method
        H = compute_homography_dlt(p1, p2)
        print("Using DLT homography computation")
    
    print("Computed Homography Matrix:")
    print(H)
    
    # Check homography quality by projecting points
    p2_projected = p2t(H, p2)
    error = np.sqrt(np.sum((p1 - p2_projected)**2, axis=0))
    print(f"Projection errors: {error}")
    print(f"Mean error: {np.mean(error):.2f} pixels")
    
    if np.mean(error) > 10:
        print("Warning: High projection error. Check point correspondences.")
    
    # Warp second image to align with first image
    img2_warped, bbox, valid_mask = imwarp_improved(img2, H, 'linear')
    
    print(f"Warped image shape: {img2_warped.shape}")
    print(f"Bounding box: {bbox}")
    
    # Show warped image for inspection
    plt.figure(figsize=(10, 6))
    plt.imshow(img2_warped.astype(np.uint8), cmap='gray')
    plt.title('Warped Image 2 (should align with Image 1)')
    plt.axis('off')
    plt.show()
    
    # Compute mosaic dimensions
    h1, w1 = img1.shape
    h2, w2 = img2_warped.shape
    
    # Calculate the overall bounding box
    x_min = min(bbox[0], 0)
    x_max = max(bbox[2], w1 - 1)
    y_min = min(bbox[1], 0)
    y_max = max(bbox[3], h1 - 1)
    
    # New mosaic dimensions
    mosaic_h = int(y_max - y_min + 1)
    mosaic_w = int(x_max - x_min + 1)
    
    print(f"Mosaic dimensions: {mosaic_h} x {mosaic_w}")
    
    # Calculate offsets for placing images
    offset_x1 = int(-x_min) if x_min < 0 else 0
    offset_y1 = int(-y_min) if y_min < 0 else 0
    offset_x2 = int(bbox[0] - x_min)
    offset_y2 = int(bbox[1] - y_min)
    
    print(f"Image 1 offset: ({offset_x1}, {offset_y1})")
    print(f"Image 2 offset: ({offset_x2}, {offset_y2})")
    
    # Create mosaic canvas
    mosaic = np.zeros((mosaic_h, mosaic_w))
    weight_map = np.zeros((mosaic_h, mosaic_w))
    
    # Place first image
    y1_start = offset_y1
    y1_end = y1_start + h1
    x1_start = offset_x1
    x1_end = x1_start + w1
    
    mosaic[y1_start:y1_end, x1_start:x1_end] = img1
    weight_map[y1_start:y1_end, x1_start:x1_end] = 1
    
    # Place second image
    y2_start = offset_y2
    y2_end = min(y2_start + h2, mosaic_h)
    x2_start = offset_x2
    x2_end = min(x2_start + w2, mosaic_w)
    
    # Ensure we don't go out of bounds
    h2_crop = y2_end - y2_start
    w2_crop = x2_end - x2_start
    
    if h2_crop > 0 and w2_crop > 0:
        # Get the valid portion of the warped image
        img2_crop = img2_warped[:h2_crop, :w2_crop]
        mask_crop = valid_mask[:h2_crop, :w2_crop]
        
        # Blend in overlapping regions
        overlap_region = weight_map[y2_start:y2_end, x2_start:x2_end] > 0
        
        # Where there's overlap, blend 50-50
        blend_mask = overlap_region & mask_crop
        mosaic[y2_start:y2_end, x2_start:x2_end] = np.where(
            blend_mask,
            0.5 * mosaic[y2_start:y2_end, x2_start:x2_end] + 0.5 * img2_crop,
            np.where(mask_crop, img2_crop, mosaic[y2_start:y2_end, x2_start:x2_end])
        )
        
        weight_map[y2_start:y2_end, x2_start:x2_end] = np.where(
            mask_crop, 1, weight_map[y2_start:y2_end, x2_start:x2_end]
        )
    
    # Display result
    plt.figure(figsize=(15, 10))
    plt.imshow(mosaic.astype(np.uint8), cmap='gray')
    plt.title('Final Mosaic Result')
    plt.axis('off')
    plt.show()

if __name__ == "__main__":
    main()

