% Analisi spettrale
% Rosa Russo VR445639
% Version that forces the refresh after each plot

clc
close all
clear all

% List of input mesh files, from .ply file
mesh_files = {'./models/Male_scale.ply', './models/Male_null.ply', './models/Male_isometric.ply'};
colors = {'g', 'b', 'r'};

figure('Name', 'Check Subplot Meshes and Spectra', 'Position', [200, 200, 1800, 400]);

vals = cell(1,3);
mesh_titles = {'Scale Shape', 'Null Shape', 'Isometric Shape'};
spec_titles = {'Scale Spectra', 'Null Spectra', 'Isometric Spectra'};

for i = 1:3
    mesh = ply_read(mesh_files{i});
    vertex = [mesh.vertex.x mesh.vertex.y mesh.vertex.z];
    triangle = zeros(length(mesh.face.vertex_indices),3);
    for j = 1:length(mesh.face.vertex_indices)
        triangle(j,:) = mesh.face.vertex_indices{j} + 1;
    end

    % Plot the mesh in the i-th subplot of the first 3 columns
    subplot(1,6,i);
    trisurf(triangle, vertex(:,1), vertex(:,2), vertex(:,3), 'FaceColor', colors{i});
    axis equal;
    axis off;
    title(mesh_titles{i});
    % Force MATLAB to render the plot
    drawnow;

    % Laplacian and autovalues
    % Compute Laplace-Beltrami operator matrices W (stiffness) and A (mass)
    % W, the Laplacian (2nd spatial derivative) of an irregular triangular mesh
    % A, the linear distances between vertices of 'face'.
    % W and A are square, [Nvertices,Nvertices] in size, sparse in nature.
    [W,A] = mesh_laplacian(vertex, triangle);
    % Compute the first 200 eigenvalues and eigenvectors of the generalized eigenproblem W*v = lambda*A*v
    [V,D] = eigs(W, A, 200, -1e-5); % The parameter sigma = -1e-5 tells eigs to find the eigenvalues closest to -1e-5, helping catch the small positive eigenvalues closest to zero
    % Store absolute values of eigenvalues for plotting the spectral signature
    vals{i} = abs(diag(D));
end

% Plot eigenvalue spectra in the next 3 subplots (columns 4 to 6)
for i = 1:3
    subplot(1,6,3+i);
    stem(vals{i}, colors{i}, 'filled', 'MarkerSize', 3); grid on;
    xlim([0 length(vals{i})]);
    title(spec_titles{i});
    drawnow;
end

