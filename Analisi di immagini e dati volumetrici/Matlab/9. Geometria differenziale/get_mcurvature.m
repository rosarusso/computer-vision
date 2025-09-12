function [m_curv] = get_mcurvature(V, T, N, L)
% Input:
%           V vertex
%           T triangles
%           N normals
%           L LB matrix
% Output: 
%           N normals
%
% 

%compute mean curvature:
m_curv = 0.5 * sum(N .* (L * V), 2);	

end

