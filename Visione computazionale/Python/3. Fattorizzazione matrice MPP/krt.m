function [K, R, t] = krt(P)

[Q, U] = qr(inv(P(1:3, 1:3)));

if U(3, 3) < 0
    U = -1 * U;
end
if U(1, 1) < 0
     E= [-1     0     0
         0      1     0
         0      0     1];
     U = E * U;
     Q = Q * E;
end
if U(2, 2) < 0
    E= [1     0     0
        0    -1     0
        0     0     1];
    U = E * U;
    Q = Q * E;
end
if det(Q) < 0 
    Q = -Q;
end

s = det(Q);
R = s * Q';
t = s * U * P(1:3, 4);
K = inv(U ./ U(3, 3));