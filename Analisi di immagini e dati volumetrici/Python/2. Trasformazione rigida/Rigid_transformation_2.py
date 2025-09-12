import open3d as o3d
import numpy as np

# Load data_i and model_i from Corr3D.mat
from scipy.io import loadmat
data = loadmat('Corr3D.mat')
data_i = data['data_i']
model_i = data['model_i']

# Convert to Open3D point clouds
source = o3d.geometry.PointCloud()
source.points = o3d.utility.Vector3dVector(data_i)

target = o3d.geometry.PointCloud()
target.points = o3d.utility.Vector3dVector(model_i)

# Subsample (optional)
source = source.random_down_sample(1/3)
target = target.random_down_sample(1/3)

# Run ICP
threshold = 0.01  # Distance threshold
trans_init = np.eye(4)  # Initial transformation

reg_result = o3d.pipelines.registration.registration_icp(
    source, target, threshold, trans_init,
    o3d.pipelines.registration.TransformationEstimationPointToPoint()
)

print("Transformation matrix:")
print(reg_result.transformation)

# Apply transformation
source.transform(reg_result.transformation)

# Visualize
o3d.visualization.draw_geometries([target, source])

