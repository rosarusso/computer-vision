import pymeshlab

# Load the registered cloud point
ms = pymeshlab.MeshSet()
ms.load_new_mesh('PCregistered.ply')

# Apply Poisson reconstruction to generate the mesh
ms.generate_surface_reconstruction_vcg()

# Save the resulting mesh
ms.save_current_mesh('mesh_completa.ply')

print("Mesh completa generata e salvata come 'mesh_completa.ply'")

