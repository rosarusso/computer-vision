# Repository for "[4S00079] Visione Computazionale" and "[4S001409] Analisi di Immagini e Dati Volumetrici" exams

This repository contains all the essential materials required to complete the exams for the courses **[4S00079] Visione Computazionale** and **[4S001409] Analisi di Immagini e Dati Volumetrici**.

Included within the repository are the assigned MATLAB scripts, which have been carefully translated into Python to facilitate broader accessibility and usability. Alongside the codebase, a literature review on 3D virtual garment modeling is provided. This review explores some of the most commonly used approaches for retrieving 3D meshes that represent garments, with the aim of identifying the most suitable method to achieve realistic clothing simulation across diverse body shapes in motion.
The report begins by presenting an overview of selected works (mostly from 2015 to 2021) that utilize 2D image data to reconstruct 3D meshes representing garments. These studies propose various innovative approaches to infer three-dimensional clothing geometry directly from images. Subsequently, the focus shifts towards methodologies that employ 3D data or volumetric inputs, which are generally more conducive to generating accurate and detailed 3D garment models.
The goal is to leverage one of these state-of-the-art mesh reconstruction techniques as a foundation for simulating realistic behavior of garments during human motion, thereby advancing the fidelity of virtual clothing simulations.


## MATLAB and Python exercises for the Course on Analisi di immagini e dati volumetrici [4S001409]
* Acquisition
  Conversion of a range map into a 3D point cloud.
  ![](./Analisi di immagini e dati volumetrici/Cloudpoint.jpg)

*Rigid Transformation
  Implementation of the Procrustes method to estimate the rigid transformation between two 3D point clouds.
  ![](<./Analisi di immagini e dati volumetrici/3D point cloud registration.jpg>)


* Pre-alignment with PCA
  Implementation of a rigid transformation estimation method that aligns the principal axes of each view using Principal Component Analysis.
  ![](<./Analisi di immagini e dati volumetrici/PCA registration.jpg>)

* Pairwise 3D Registration
  Implementation of the 3D registration algorithm for a pair of point clouds using the Iterative Closest Point (ICP) algorithm.
  ![](<./Analisi di immagini e dati volumetrici/3D point cloud registration.jpg>)

* Registration with X-84 ICP Algorithm
  Implementation of the robust variant of the ICP algorithm where closest points are filtered by automatically discarding outliers using the X-84 rule.
  ![](<./Analisi di immagini e dati volumetrici/ICP alignment comparison - py_mat.jpg>)

* Multiview Registration
  Implementation of the basic method for global registration by combining pairwise registrations.
  ![](<./Analisi di immagini e dati volumetrici/Registered mesh.png>)

* Polygon Mesh Generation from a Range Map
  Define a triangular mesh starting from the point cloud obtained from a range map.
  ![](<./Analisi di immagini e dati volumetrici/Depth image - 3D point cloud.jpg>)

* Complete Object Mesh Generation
  Define a single polygonal mesh starting from a set of 3D point clouds represented by aligned partial views.
  ![](<./Analisi di immagini e dati volumetrici/AlphaShape Reconstruction with Rotation.png>)

* Differential Geometry
  Estimate the mean curvature from the combination of vertex normals and the Laplacian operator associated with the coordinate function.
  ![](<./Analisi di immagini e dati volumetrici/Original mesh and mean curvature-py.jpg>)

* Spectral Analysis of 3D Shapes
  Calculate eigenvalues and eigenfunctions of the Laplace–Beltrami Operator and plot the sequence of eigenvalues (i.e., ShapeDNA).
  ![](<./Analisi di immagini e dati volumetrici/Meshes and Eigenvalue spectra.png>)