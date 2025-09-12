function writeply(verts, faces, filename, mode)
%WRITEPLY Write data to PLY file.
%
%   writeply(verts, faces, filename, mode) writes the triangular mesh
%   specified by verts and faces into the PLY file filename. The mode
%   attribute can be 'ascii', 'binary_little_endian' or 'binary_big_endian'
%   (DEFAULT: 'binary_big_endian')

% Author: G. Peyr, 2004
% Edit: M. Schivi, 2015

if nargin < 4
    mode = 'binary_big_endian';
end

% Check inputs
assert(ismatrix(verts) && size(verts, 2) == 3, ...
    'WRITEPLY:InvalidInput', 'verts must be a N x 3 matrix');

assert(ismatrix(faces) && size(faces, 2) == 3, ...
    'WRITEPLY:InvalidInput', 'faces must be a M x 3 matrix');

% Add verts to data
data.vertex.x = verts(:,1);
data.vertex.y = verts(:,2);
data.vertex.z = verts(:,3);

% Add faces to data
data.face.vertex_indices = mat2cell(faces - 1, ...
    ones(size(faces, 1), 1), size(faces, 2));

plywrite(data, filename, mode);
