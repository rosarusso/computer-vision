function R = eul(a)

% EUL: compute a rotation matrix given euler angles (yaw,pitch,roll)

phi   = a(3);
theta = a(2);
psi   = a(1);

R = [
    cos(phi)  -sin(phi) 0
    sin(phi)   cos(phi) 0
    0             0     1]*[
    cos(theta)  0 sin(theta)
    0         1    0
    -sin(theta) 0 cos(theta)]*[
    1         0    0
    0  cos(psi)  -sin(psi)
    0  sin(psi)   cos(psi)];
