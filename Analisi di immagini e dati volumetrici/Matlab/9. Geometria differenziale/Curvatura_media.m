% Stima della curvatura media in ogni vertice data una mesh poligonale in
% input, ricavata dalla combinazione tra la matrice di Laplace-Beltrami e
% le normali di ogni vertice
% Rosa Russo VR445639

clc
close all
clear all

% Getting the mesh from .ply file
mesh = ply_read('./models/Centaur_null.ply');

% Vertices
vertex=[mesh.vertex.x mesh.vertex.y mesh.vertex.z];
% Extract face triangles (we add 1 because MATLAB uses 1-based indexing)
triangle=zeros(length(mesh.face.vertex_indices),3);
for i=1:length(mesh.face.vertex_indices)
    triangle(i,:)=mesh.face.vertex_indices{i}+1;
end

figure;
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), ones(length(vertex), 1)), shading interp, axis image, axis tight, lighting phong, camlight

% Compute Laplace-Beltrami matrix W and mass matrix A
[W,A]=mesh_laplacian(vertex,triangle);

L=inv(A)*W; % Laplacian operator L=A^(-1)*W

% Compute normals at each vertex
N=getNormals(vertex,triangle);
% Estimate mean curvature at each vertex
curv=get_mcurvature(vertex,triangle,N,L);

figure;
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), curv), shading interp, axis image, axis tight, lighting phong, camlight
colorbar;
%{
The colorbar shows the color scale used to represent mean curvature values on the mesh vertices.
Colors towards the cooler end of the scale indicate low curvature regions,
while warmer colors indicate high curvature regions.
%}
title('Estimated mean curvature');

colorbar.Label.String = 'Mean Curvature Value';
colorbar.Label.FontSize = 12;
colorbar.Label.FontWeight = 'bold';
