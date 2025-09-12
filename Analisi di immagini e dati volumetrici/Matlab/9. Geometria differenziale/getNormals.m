function N = getNormals(V, T)
% Input:
%           V vertex
%           T triangles
% Output: 
%           N normals
%
%
	N = cross( V(T(:,1),:) - V(T(:,2),:), V(T(:,1),:) - V(T(:,3),:));
	N = [accumarray(T(:), repmat(N(:,1), [3,1])) , accumarray(T(:), repmat(N(:,2) , [3,1])), accumarray(T(:), repmat(N(:,3), [3,1]))];
	N = N ./ repmat(sqrt(sum(N.^2, 2)), [1, 3]);
end