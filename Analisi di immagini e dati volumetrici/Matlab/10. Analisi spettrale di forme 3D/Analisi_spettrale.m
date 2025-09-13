% Analisi spettrale
% Rosa Russo VR445639

clc
close all
clear all

% First shape
% Autovalues and autovectors extraction

% First mesh loading from .ply file
mesh = ply_read('./models/Male_scale.ply');

% Vertices
vertex = [mesh.vertex.x mesh.vertex.y mesh.vertex.z];
% Triangles
triangle = zeros(length(mesh.face.vertex_indices),3);
for i=1:length(mesh.face.vertex_indices)
    triangle(i,:) = mesh.face.vertex_indices{i}+1;
end

% Figure and subplot setup
figure('Name','Spectral Analysis', 'Position', [200, 200, 1200, 600]);
subplot(2,3,1);
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), ones(length(vertex), 1)), shading interp, axis image, axis tight, lighting phong, camlight
title('Scale shape');

% W, the Laplacian (2nd spatial derivative) of an irregular triangular mesh
% A, the linear distances between vertices of 'face'.
% W and A are square, [Nvertices,Nvertices] in size, sparse in nature.
[W,A] = mesh_laplacian(vertex,triangle);

% Autovectors and autovalues decomposition
[vet val] = eigs(W,A,200,-1e-5);
vals1 = abs(diag(val));

% Spectrum
%figure(2);
%plot(abs(val), 'gx');
%grid on;
%hold on;
%ylabel('eigenvalues');

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

subplot(2,3,2);
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), ones(length(vertex), 1)), shading interp, axis image, axis tight, lighting phong, camlight
title('Null shape');

% W, the Laplacian (2nd spatial derivative) of an irregular triangular mesh
% A, the linear distances between vertices of 'face'.
% W and A are square, [Nvertices,Nvertices] in size, sparse in nature.
[W,A] = mesh_laplacian(vertex,triangle);

% Autovectors and autovalues decomposition
[vet val]=eigs(W,A,200,-1e-5);
vals2=abs(diag(val));

% Spectrum
%figure(4);
%plot(abs(val), 'bo');


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

subplot(2,3,3);
trisurf(triangle,vertex(:,1), vertex(:,2), vertex(:,3), ones(length(vertex), 1)), shading interp, axis image, axis tight, lighting phong, camlight
title('Isometric shape');

% Laplacian and autovalues
% W, the Laplacian (2nd spatial derivative) of an irregular triangular mesh
% A, the linear distances between vertices of 'face'.
% W and A are square, [Nvertices,Nvertices] in size, sparse in nature.
[W,A] = mesh_laplacian(vertex,triangle);

% Autovectors and autovalues decomposition
[vet val]=eigs(W,A,200,-1e-5);
vals3=abs(diag(val));

% Spectrum
%figure(6);
%plot(abs(val), 'r+');
%legend('scale','null','isometric');

% Spectrum, all of the three
subplot(2,3,4);
stem(vals1, 'g', 'filled', 'MarkerSize', 3);
title('Scale Spectra');
grid on;

subplot(2,3,5);
stem(vals2, 'b', 'filled', 'MarkerSize', 3);
title('Null Spectra');
grid on;

subplot(2,3,6);
stem(vals3, 'r', 'filled', 'MarkerSize', 3);
title('Isometric Spectra');
grid on;

