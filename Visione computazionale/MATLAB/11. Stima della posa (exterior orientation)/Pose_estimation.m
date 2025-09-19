clear all
close all
clc

% Camera information loading
load('imgInfo.mat')

img = imread('cav.jpg');

% Image information
p2D = imgInfo.punti2DImg;
p3D = imgInfo.punti3DImg;
p3Dn = [p3D'; ones(1, size(p3D, 1))]; % homogeneous
K = imgInfo.K;

N = 100; % Number of points to be used

% 3D points
figure(1)
scatter3(p3D(:,1),p3D(:,2),p3D(:,3),5,'c');
axis equal
% 3D points projection (with calibration info)
figure(2)
imshow(img);
hold on;
plot(p2D(:,1), p2D(:,2),'r.');
P = K*[imgInfo.R imgInfo.T];
[u,v] = proj(P,p3D);
plot(u,v,'go');
ps = normalize(P*p3Dn);
plot(ps(1,:),ps(2,:),'go');
% 3D points
figure(3)
scatter3(p3D(:,1),p3D(:,2),p3D(:,3),5,'c');
axis equal

% Exterior orientation
[G,~] = exterior_fiore(K, p3D(1:N, :)', p2D(1:N, :)');

% 3D points projection with Fiore externals extimation
figure(4);
imshow(img);
hold on;
plot(p2D(:,1), p2D(:,2),'r.');
P1 = K*G;
[u1,v1] = proj(P1,p3D);
plot(u1,v1,'bo');
ps1 = normalize(P1*p3Dn);
plot(ps1(1,:),ps1(2,:),'bo');

% Comparison between the pose estimated using the calibration and the one
% estimated using the Fiore method
[imgInfo.R imgInfo.T G]