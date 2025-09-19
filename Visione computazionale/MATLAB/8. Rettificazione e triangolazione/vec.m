function v = vec(A);

% VEC puts the elements of a matrix in a column vector.
%
% v = vec(A) given a matrix A (n x m), it returns a column vector
% containing n x m elements, which are the union of the matrix A columns.

v = reshape(A,size(A,1)*size(A,2),1);
