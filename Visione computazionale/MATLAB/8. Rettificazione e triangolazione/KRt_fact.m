%factorize P into K, R and t
function [K,R,t] = KRt_fact(P)
    Q = P(:,1:3);
    q = P(:,4);
    [O,T] = qr(inv(Q));
    %internal
    K = inv(T);
    detO = det(O);
    %sign correction
    R = detO*O.';
    %translation
    t = T*(detO*q);
    %normalization
    K = K/K(3,3);
end
