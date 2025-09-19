function M = intersect_base(PPM,m)
numP = length(m{1});
numV = length(m);

M = [];
for i = 1:numP  
    A=[];
    for view = 1:numV
        PPM{view} = PPM{view}./norm(PPM{view}(3,1:3));
        A = [A
             PPM{view}(1,:) - m{view}{i}(1)*PPM{view}(3,:)
             PPM{view}(2,:) - m{view}{i}(2)*PPM{view}(3,:)];
    end
    x = ns(A);
    M = [M x(1:3)./x(4)];
end

