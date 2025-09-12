function plywrite(elements, path, format, prec)
%PLYWRITE Write 3D data as a PLY file.
%   plywrite(elements, path) writes the structure elements as a binary PLY 
%   file. Every field of elements is interpreted as an element and every 
%   subfield as an element property. Each subfield of property data must 
%   either be an array or a cell array of arrays. All property data in an 
%   element must have the same length.
%
%   A common PLY data structure has the following fields:
%       elements.vertex.x = x coordinates, [Nx1] real array
%       elements.vertex.y = y coordinates, [Nx1] real array
%       elements.vertex.z = z coordinates, [Nx1] real array
%
%       elements.face.vertex_indices = vertex index lists,
%           an {Mx1} cell array where each cell holds a one-dimesional 
%           array (of any length) of vertex indices.
%
%   Some other common data fields:
%       elements.vertex.nx = x coordinate of normal, [Nx1] real array
%       elements.vertex.ny = y coordinate of normal, [Nx1] real array
%       elements.vertex.nz = z coordinate of normal, [Nx1] real array
%
%       elements.edge.vertex1 = index to a vertex, [Px1] integer array
%       elements.edge.vertex2 = second vertex index, [Px1] integer array
%
%   Many other fields and properties can be added. The PLY format is not 
%   limited to the naming in the examples above -- they are only the 
%   conventional naming.
%
%   plywrite(..., format) write the PLY with a specified data format, where 
%   format is:
%       'ascii'                 ASCII text data
%       'binary_little_endian'  binary data, little endian
%       'binary_big_endian'     binary data, big endian (default)
%
%   plywrite(..., 'double') write floating-point data as double precision 
%   rather than in the default single precision.
%
%   See also PLYREAD.

% Pascal Getreuer 2004
% Edit: M. Schivi 2015

% Check the inputs
if nargin < 4
    % No precision, use single
    prec = 'single';
    
    if nargin < 3
        % No format and precision
        format = 'binary_big_endian';
    elseif strcmpi(format, 'double')
        prec = 'double';
        format = 'binary_big_endian';
    end
end

% Open the output file
[fid, msg] = fopen(path, 'wt');

if fid == -1
    error(msg);
end

plyTypes        = {...
    'int', 'float', 'double', ...                   % Array types
    'int', 'float', 'double'};                      % Cell types
fwriteTypes     = {'int32', 'single', 'double'};
matlabTypes     = {'int32', 'single', 'double'};
printfFormat    = {'%d', '%-.6f', '%-.14e'};

% Write PLY header
fprintf(fid, 'ply\n');
fprintf(fid, 'format %s 1.0\n', format);

% Retrieve element names
elementNames = fieldnames(elements);
nElements = length(elementNames);
data = cell(nElements, 1);

% For each element
for i = 1:nElements
    % Is struct?
    tmp = eval(['isstruct(elements.', elementNames{i}, ');']);
    
    if tmp
        propertyNames{i} = eval(...
            ['fieldnames(elements.', elementNames{i}, ');']);
    else
        propertyNames{i} = [];
    end
    
    if ~isempty(propertyNames{i})
        % Store data to write
        data{i}{1} = eval(...
            ['elements.', elementNames{i}, '.', propertyNames{i}{1}, ';']);
        elementCount(i) = numel(data{i}{1});
        type{i} = zeros(length(propertyNames{i}), 1);
    else
        elementCount(i) = 0;
    end
    
    fprintf(fid, 'element %s %u\n', elementNames{i}, elementCount(i));
    
    % For each property
    for j = 1:length(propertyNames{i})
        % Store data to write
        data{i}{j} = eval(...
            ['elements.', elementNames{i}, '.', propertyNames{i}{j}, ';']);
        
        if elementCount(i) ~= numel(data{i}{j});
            fclose(fid);
            error('All property data in element must have the same length');
        end
        
        if iscell(data{i}{j})
            type{i}(j) = 3;
            data{i}{j} = data{i}{j}{1};
        end
        
        for k = 1:length(matlabTypes)
            if isa(data{i}{j}, matlabTypes{k})
                type{i}(j) = type{i}(j) + k;
                break;
            end
        end
        
        % Try to convert float data to integer data
        if type{i}(j) <= 3      % Array data
            if any(strcmp({'single', 'double'}, matlabTypes{type{i}(j)}))
                if ~any(floor(data{i}{j}) ~= data{i}{j})        
                    % data is integer
                    type{i}(j) = 1;
                end
            end
        else                    % Cell array data
            data{i}{j} = eval(...
                ['elements.', elementNames{i}, '.', propertyNames{i}{j}, ';']);
            tmp = 1;
            
            for k = 1:numel(data{i}{j})
                tmp = tmp & all(floor(data{i}{j}{k}) == data{i}{j}{k});
            end
            
            type{i}(j) = tmp + 3;
        end
    
        % Convert double to single if specified
        if rem(type{i}(j), 3) == 0 && strcmpi(prec, 'single')
            type{i}(j) = type{i}(j) - 1;
        end
    
        if type{i}(j) <= 3
            fprintf(fid, 'property %s %s\n', ...
                plyTypes{type{i}(j)},...
                propertyNames{i}{j});
        else
            fprintf(fid, 'property list uchar %s %s\n', ...
                plyTypes{type{i}(j) - 3}, ...
                propertyNames{i}{j});
        end
    end
end

fprintf(fid, 'end_header\n');

switch format
    case 'ascii'
        format = 0;
    case 'binary_little_endian'
        fclose(fid);
        fid = fopen(path, 'a', 'ieee-le');
        format = 1;
    case 'binary_big_endian'
        fclose(fid);
        fid = fopen(path, 'a', 'ieee-be');
        format = 2;
end

for i = 1:nElements
    if ~isempty(propertyNames{i})
        if ~format      % Write ASCII data
            for k = 1:elementCount(i)
                for j = 1:length(propertyNames{i})
                    if type{i}(j) < 3
                        fprintf(fid, ...
                            [printfTypeChar{Type{i}(j)}, ' '], ...
                            data{i}{j}(k));
                    else
                        fprintf(fid, ...
                            '%u%s ', ...
                            length(data{i}{j}{k}), ...
                            sprintf(...
                            [' ', printfTypeChar{type{i}(j) - 3}], ...
                            data{i}{j}{k}));
                    end
                end
            
                fprintf(fid, '\n');
            end
        else            % Write binary data
            if all(type{i} < 3) && all(type{i} == type{i}(1))
                % Property data without list types (fast)
                tmp = zeros(length(propertyNames{i}), elementCount(i));
         
                for j = 1:length(propertyNames{i})
                    tmp(j, :) = data{i}{j}(:)';
                end
         
                fwrite(fid, tmp, fwriteTypes{type{i}(j)});
            elseif all(type{i} >= 3)
                % Only list types
                type{i} = type{i} - 3;
            
                if length(propertyNames{i}) == 1
                    % Only one list property
                    tmp = fwriteTypes{type{i}(1)};
            
                    for k = 1:elementCount(i)
                        fwrite(fid, length(data{i}{1}{k}), 'uchar');
                        fwrite(fid, data{i}{1}{k}, tmp);
                    end
                else
                % Multiple list properties
                    for k = 1:elementCount(i)
                        for j = 1:length(propertyNames{i})
                            fwrite(fid, ...
                                length(data{i}{j}{k}), 'uchar');
                            fwrite(fid, ...
                                data{i}{j}{k}, fwriteTypes{type{i}(j)});
                        end
                    end
                end
            else
                % Mixed type
                for k = 1:elementCount(i)
                    for j = 1:length(propertyNames{i})
                        if type{i}(j) < 3
                            fwrite(fid, ...
                                data{i}{j}(k), fwriteTypes{type{i}(j)});
                        else
                            fwrite(fid, ...
                                length(data{i}{j}{k}), 'uchar');
                            fwrite(fid, ...
                                data{i}{j}{k}, fwriteTypes{Type{i}(j) - 3});
                        end
                    end
                end
            end
        end
    end
end

fclose(fid);
