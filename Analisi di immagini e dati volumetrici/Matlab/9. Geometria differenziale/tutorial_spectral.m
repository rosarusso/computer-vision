% 
% Eurographics 2012 - Tutorial on "Diffusion Geometry in Shape Analysis"
%
%  Umberto Castellani, Alex Bronstein, Michael Bronstein
%

clc
close all
clear

% Read model as ply format
%
data = ply_read('./models/Male_null.ply');


vertex=[data.vertex.x data.vertex.y data.vertex.z];
%triangoli
triangle=zeros(length(data.face.vertex_indices),3);
for i=1:length(data.face.vertex_indices)
    triangle(i,:)=data.face.vertex_indices{i}+1;
end
%
shape.X = data.vertex.x;
shape.Y = data.vertex.y;
shape.Z = data.vertex.z;
TRIV = zeros(length(data.face.vertex_indices),3);
for k = 1:length(data.face.vertex_indices)
   TRIV(k,:) = data.face.vertex_indices{k} + 1; 
end
shape.TRIV = TRIV;

% Show 3D model
%
figure(1)
trisurf(shape.TRIV,shape.X,shape.Y,shape.Z, ones(length(shape.X), 1)), shading interp, axis image, axis tight, lighting phong, camlight


% Eigendecomposition
%
% Compute matrices (see http://www.geomtop.org/software/meshlp.html)
%
%[W, A] = mshlp_matrix(shape);
%Am = sparse([1:length(A)], [1:length(A)], A);
[W,Am]=mesh_laplacian(vertex,triangle);
warning off;
% Compute eigendecomposition
neig=200;
[evecs evals] = eigs(W, Am, neig, -1e-5);
warning on;
evals = abs(diag(evals));

% Show eigenfunction
viseig=1;

if(viseig)
    for i=2:20:100
        figure
        trisurf(shape.TRIV,shape.X,shape.Y,shape.Z,evecs(:,i)), shading interp, axis image, axis tight, lighting phong, camlight
        colorbar
    end
end

input('...');
% show spectra:
figure(101)
plot(abs(evals), 'rx');
grid on
ylabel('eigenvalues');
disp('This is the spectra of a basic shape');
input('...');



% load other model - the same model with a isometic variation
data = ply_read('./models/Male_isometric.ply');
shape.X = data.vertex.x;
shape.Y = data.vertex.y;
shape.Z = data.vertex.z;
TRIV = zeros(length(data.face.vertex_indices),3);
for k = 1:length(data.face.vertex_indices)
   TRIV(k,:) = data.face.vertex_indices{k} + 1; 
end
shape.TRIV = TRIV;

vertex=[data.vertex.x data.vertex.y data.vertex.z];
%triangoli
triangle=zeros(length(data.face.vertex_indices),3);
for i=1:length(data.face.vertex_indices)
    triangle(i,:)=data.face.vertex_indices{i}+1;
end

% Show 3D model
%
figure(200)
trisurf(shape.TRIV,shape.X,shape.Y,shape.Z, ones(length(shape.X), 1)), shading interp, axis image, axis tight, lighting phong, camlight
%[W, A] = mshlp_matrix(shape);
%Am = sparse([1:length(A)], [1:length(A)], A);
[W,Am]=mesh_laplacian(vertex,triangle);

warning off;
% Compute eigendecomposition
neig=200;
[evecs_iso evals_iso] = eigs(W, Am, neig, -1e-5);
warning on;
evals_iso = abs(diag(evals_iso));
%
% Show spectra
figure(101)
hold on
plot(abs(evals_iso), 'bo');
ylabel('eigenvalues');
grid on
disp('This is a spectra of the same model with an isometric transformation');
disp('The spectra is maintained.');
input('...');


% load other model
data = ply_read('./models/Male_scale.ply');
shape.X = data.vertex.x;
shape.Y = data.vertex.y;
shape.Z = data.vertex.z;
TRIV = zeros(length(data.face.vertex_indices),3);
for k = 1:length(data.face.vertex_indices)
   TRIV(k,:) = data.face.vertex_indices{k} + 1; 
end
shape.TRIV = TRIV;

vertex=[data.vertex.x data.vertex.y data.vertex.z];
%triangoli
triangle=zeros(length(data.face.vertex_indices),3);
for i=1:length(data.face.vertex_indices)
    triangle(i,:)=data.face.vertex_indices{i}+1;
end

% Show 3D model
%
figure(300)
trisurf(shape.TRIV,shape.X,shape.Y,shape.Z, ones(length(shape.X), 1)), shading interp, axis image, axis tight, lighting phong, camlight
%[W, A] = mshlp_matrix(shape);
%Am = sparse([1:length(A)], [1:length(A)], A);
[W,Am]=mesh_laplacian(vertex,triangle);
warning off;
% Compute eigendecomposition
neig=200;
[evecs_scale evals_scale] = eigs(W, Am, neig, -1e-5);
warning on;
evals_scale = abs(diag(evals_scale));
%
% Show spectra
figure(101)
hold on
plot(abs(evals_scale), 'm+');
grid on
ylabel('eigenvalues');
legend('null-shape', 'isometric', 'scale');
disp('Now the spectra of the same model with a scale variation is added');
disp('The spectra is different (i.e., no scale invariance).');



