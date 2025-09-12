function mesh = readply(filename)
%READPLY Read data from PLY file.
%
%   mesh = readply(filename) reads the triangular mesh specified
%   by the PLY file filename and retuns it.
%
%   See also TRIANGULARMESH.

% Author: G. Peyr, 2003
% Edit: M. Schivi, 2015 

% Read file
data = plyread(filename);

% Check fields
assert(isfield(data, 'vertex') && isfield(data, 'face'), ...
    'READPLY:InvalidInput', 'ply file must have fields vertex and face.')

% Check triangular
lengths = cellfun('length', data.face.vertex_indices);
assert(all(lengths == 3), ...
    'READPLY:InvalidInput', 'Mesh must be triangular.');

% Get verts
verts = [data.vertex.x, data.vertex.y, data.vertex.z];

% Get faces
faces = cell2mat(data.face.vertex_indices) + 1;

% Get TriangularMesh
mesh = TriangularMesh(verts, faces);
