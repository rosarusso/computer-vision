% Analisi spettrale
% Rosa Russo VR445639

clc
close all
clear all

% First shape
% Autovalues and autovectors extraction

% Mesh loading from .ply file
mesh = ply_read('./models/Male_scale.ply');

% Vertices
vertex = [mesh.vertex.x mesh.vertex.y mesh.vertex.z];
% Triangles
triangle = zeros(length(mesh.face.vertex_indices),3);
for i=1:length(mesh.face.vertex_indices)
    triangle(i,:) = mesh.face.vertex_indices{i}+1;
end

figure(2);
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), ones(length(vertex), 1)), shading interp, axis image, axis tight, lighting phong, camlight

% W, the Laplacian (2nd spatial derivative) of an irregular triangular mesh
% A, the linear distances between vertices of 'face'.
% W and A are square, [Nvertices,Nvertices] in size, sparse in nature.
[W,A] = mesh_laplacian(vertex,triangle);

% Autovectors and autovalues decomposition
[vet val] = eigs(W,A,200,-1e-5);
val = abs(diag(val));

% Spectrum
figure(1);
plot(abs(val), 'gx');
grid on;
hold on;
ylabel('eigenvalues');

%% Second shape
% Autovalues and autovectors extraction

% Mesh loading from .ply file
mesh = ply_read('./models/Male_null.ply');

% Vertices
vertex = [mesh.vertex.x mesh.vertex.y mesh.vertex.z];
% Triangles
triangle = zeros(length(mesh.face.vertex_indices),3);
for i=1:length(mesh.face.vertex_indices)
    triangle(i,:) = mesh.face.vertex_indices{i}+1;
end

figure;
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), ones(length(vertex), 1)), shading interp, axis image, axis tight, lighting phong, camlight

% W, the Laplacian (2nd spatial derivative) of an irregular triangular mesh
% A, the linear distances between vertices of 'face'.
% W and A are square, [Nvertices,Nvertices] in size, sparse in nature.
[W,A] = mesh_laplacian(vertex,triangle);

% Autovectors and autovalues decomposition
[vet val]=eigs(W,A,200,-1e-5);
val=abs(diag(val));

% Spectrum
figure(1);
plot(abs(val), 'bo');


%% Third shape
% Autovalues and autovectors extraction

% Mesh loading from .ply file
mesh = ply_read('./models/Male_isometric.ply');

% Vertices
vertex = [mesh.vertex.x mesh.vertex.y mesh.vertex.z];
% Triangles
triangle = zeros(length(mesh.face.vertex_indices),3);
for i=1:length(mesh.face.vertex_indices)
    triangle(i,:) = mesh.face.vertex_indices{i}+1;
end

figure;
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), ones(length(vertex), 1)), shading interp, axis image, axis tight, lighting phong, camlight

% W, the Laplacian (2nd spatial derivative) of an irregular triangular mesh
% A, the linear distances between vertices of 'face'.
% W and A are square, [Nvertices,Nvertices] in size, sparse in nature.
[W,A] = mesh_laplacian(vertex,triangle);

% Autovectors and autovalues decomposition
[vet val]=eigs(W,A,200,-1e-5);
val=abs(diag(val));

% Spectru,
figure(1);
plot(abs(val), 'r+');
legend('scale','null','isometric');