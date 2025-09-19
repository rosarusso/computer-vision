%nullspace
function v = ns(A)
    [~, ~, V] = svd(A);
    v = V(:,end);
end


