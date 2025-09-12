% Curvature driver
%
% Umberto castellani 2019

clc
close all
clear

% Read model as ply format
%
data = ply_read('./models/Centaur_null.ply');
%Vertices
V=[data.vertex.x, data.vertex.y, data.vertex.z];
%Triangles
TRIV = zeros(length(data.face.vertex_indices),3);
for k = 1:length(data.face.vertex_indices)
   TRIV(k,:) = data.face.vertex_indices{k} + 1;
end


% Show 3D model
%
figure(1)
trisurf(TRIV,V(:,1), V(:,2), V(:,3), ones(length(V), 1)), shading interp, axis image, axis tight, lighting phong, camlight

% Prepare Shape structure:
shape.X = data.vertex.x;
shape.Y = data.vertex.y;
shape.Z = data.vertex.z;
shape.TRIV=TRIV;
% Cot-weight matrix and Area matrix:
[W, A] = mshlp_matrix(shape);
Am = sparse([1:length(A)], [1:length(A)], A);
% LB-matrix
L=inv(Am)*W;

% Compute normals
N = getNormals(V, TRIV);

%Compute mean curvature
[m_curv] = get_mcurvature(V, TRIV, N, L);
%
figure(2)
trisurf(TRIV,V(:,1), V(:,2), V(:,3),m_curv), shading interp, axis image, axis tight, lighting phong, camlight
colorbar

