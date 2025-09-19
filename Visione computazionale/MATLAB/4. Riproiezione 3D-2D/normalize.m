%normalize points in homogeneous coordinates
function m = normalize(m)
    s = size(m);
    m = m./m(s(1),:);
end