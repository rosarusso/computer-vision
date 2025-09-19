%fundamental matrix from P
function [F,el,er] = fundamental(pl,pr)

% Optical centers
cl = -inv(pl(:,1:3))*pl(:,4);
cr = -inv(pr(:,1:3))*pr(:,4);

% Epipolar
%el = pl*[cr' 1]';
%er = pr*[cl' 1]';
el = pl*[cr; 1];
er = pr*[cl; 1];

el = el./norm(el);
er = er./norm(er);

% Fundamental Matrix
F = [   0    -er(3)  er(2)
     er(3)    0   -er(1)
    -er(2)  er(1)   0]*pr(:,1:3)*inv(pl(:,1:3));

F = F./norm(F);
end
