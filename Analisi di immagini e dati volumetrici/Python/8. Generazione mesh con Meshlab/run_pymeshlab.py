import pymeshlab

# Carica la nuvola di punti registrata
ms = pymeshlab.MeshSet()
ms.load_new_mesh('PCregistered.ply')

# Applica la ricostruzione di Poisson per generare la mesh
ms.generate_surface_reconstruction_vcg()

# Salva la mesh risultante
ms.save_current_mesh('mesh_completa.ply')

print("Mesh completa generata e salvata come 'mesh_completa.ply'")

