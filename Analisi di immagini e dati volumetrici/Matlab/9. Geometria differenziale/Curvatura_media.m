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
% Triangles
triangle=zeros(length(mesh.face.vertex_indices),3);
for i=1:length(mesh.face.vertex_indices)
    triangle(i,:)=mesh.face.vertex_indices{i}+1;
end

figure;
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), ones(length(vertex), 1)), shading interp, axis image, axis tight, lighting phong, camlight

[W,A]=mesh_laplacian(vertex,triangle);

L=inv(A)*W;

N=getNormals(vertex,triangle);
curv=get_mcurvature(vertex,triangle,N,L);

figure;
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), curv), shading interp, axis image, axis tight, lighting phong, camlight
colorbar;
