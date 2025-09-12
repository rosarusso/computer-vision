% Pre-allineamento con PCA
% Rosa Russo VR445639

clc
clear all
close all

% Point cloud
data = load('a4000001.cnn');

% Rotation vector
rotv = rand(3,1);

% 3D rotation matrix
Rz = [cos(rotv(3)) -sin(rotv(3)) 0;
      sin(rotv(3))  cos(rotv(3)) 0;
         0              0        1];
   
Ry = [cos(rotv(2)) 0 sin(rotv(2));
         0         1       0;
     -sin(rotv(2)) 0 cos(rotv(2))];

Rx = [1      0             0;
      0 cos(rotv(1)) -sin(rotv(1));
      0 sin(rotv(1))  cos(rotv(1))];

R = Rz*Ry*Rx;

% Rototranslation matrix
G = [R,  zeros(3,1);];

datamod = [data'; ones(size(data,1),1)'];

% Rotated data with random Gaussian noise
datarot = (G*datamod)';
datarot = datarot+20.*rand(size(datarot,1),3);

figure;
plot3(data(:,1), data(:,2), data(:,3), '.r');
hold on;
plot3(datarot(:,1), datarot(:,2), datarot(:,3), '.b');
hold on;
grid on;
axis equal;

% Clouds centroids (initial point cloud and rotated one)
C1 = mean(data);
plot3(C1(1),C1(2),C1(3),'.g','MarkerSize',40);
C2 = mean(datarot);
plot3(C2(1),C2(2),C2(3),'.g','MarkerSize',40);

% PCA computation and normalization - first view
R1 = pca(data);

for i=1:size(R1,2)
    R1(:,i) = R1(:,i)/norm(R1(:,i));
end

if (dot(cross(R1(:,1), R1(:,2)), R1(:,3))>0)
    disp('Righ-hand, OK');
else
    disp('Left-hand, righ-hand conversion');
    R1(:,3) = -R1(:,3);
end

G1 = [R1' -(R1'*C1')];

datacan = (G1*datamod)';
plot3(datacan(:,1), datacan(:,2), datacan(:,3), '.y');

% PCA computation and normalization - second view
R2 = pca(datarot);

for i=1:size(R2,2)
    R2(:,i) = R2(:,i)/norm(R2(:,i));
end

if (dot(cross(R2(:,1), R2(:,2)), R2(:,3))>0)
    disp('Righ-hand, OK');
else
    disp('Left-hand, righ-hand conversion');
    R2(:,3) = -R2(:,3);
end

% 4 possible cases checking (only right-hand cases)
R2_all = zeros(3,3,4);
R2_all(:,:,1) = R2; %(+,+,+)
R2_all(:,:,2) = [R2(:,1), -R2(:,2), -R2(:,3)]; % (+,-,-)
R2_all(:,:,3) = [-R2(:,1), R2(:,2), -R2(:,3)]; % (-,+,-)
R2_all(:,:,4) = [-R2(:,1), -R2(:,2), R2(:,3)]; % (-,-,+)

% Reference system on first view
GReftot = zeros(4,4,4);
Tr = zeros(1,4);

for i=1:4 % 4 possibilities
    R2 = R2_all(:,:,i);
    G2 = [R2' -(R2'*C2')];
    datarotmod = [datarot'; ones(size(datarot,1),1)'];
    datarot_can = (G2*datarotmod)';
    C1_can = mean(datacan);
    C2_can = mean(datarot_can);
    
    figure;
    disp(['Evaluation of case #' num2str(i)]);
    % Canonical pose:
    plot3(datacan(:,1), datacan(:,2), datacan(:,3), '.r');
    hold on
    plot3(datarot_can(:,1), datarot_can(:,2), datarot_can(:,3), '.b');
    hold on
    
    plot3(C1_can(1), C1_can(2), C1_can(3), '.r', 'Markersize', 30);
    plot3(C2_can(1), C2_can(2), C2_can(3), '.r', 'Markersize', 30);
    grid on
    axis equal
    
    figure;
    G1mod = [G1; 0 0 0 1];
    G2mod = [G2; 0 0 0 1];
    GRef = G1mod\G2mod;
    GReftot(:,:,i) = GRef;
    Tr(i) = trace(GRef(1:3,1:3));
    
    disp(trace(GRef(1:3,1:3)));
    disp(rad2deg(ieul(GRef(1:3,1:3))));
    
    % GRef can be used as pre-alignment of ICP
    datarot_ref = (GRef*datarotmod)';
    plot3(data(:,1), data(:,2), data(:,3), '.r');
    hold on
    plot3(datarot_ref(:,1), datarot_ref(:,2), datarot_ref(:,3), '.b');
    grid on
    axis equal
end

% Euristic for automatic selection (max Trace --> min rotation):
[drop,optimal] = max(Tr);
disp(['Best case according to heuristic on trace: ' num2str(optimal)])