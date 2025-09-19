%fundamental matrix from P
function [F,el,er] = fundamental(pl,pr)

%center
cl = -inv(pl(:,1:3))*pl(:,4);
cr = -inv(pr(:,1:3))*pr(:,4);

%epipolar
el = pl*[cr' 1]';
er = pr*[cl' 1]';

%fundamental
F=[   0    -er(3)  er(2)
     er(3)    0   -er(1)
    -er(2)  er(1)   0   ]*pr(:,1:3)*inv(pl(:,1:3));
end

