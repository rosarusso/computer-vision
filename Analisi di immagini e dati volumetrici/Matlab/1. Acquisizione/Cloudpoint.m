% Acquisition
% Rosa Russo VR445639
% From range data map to 3D point cloud

clc
close all
clear all

% Intrinsic matrix
DeptMat = [575.8, 0.0, 319.5;
          0.0, 575.8, 239.5;
          0.0, 0.0, 1.0];

% Coordinates
cenX = DeptMat(1,3);
cenY = DeptMat(2,3);
focX = DeptMat(1,1);
focY = DeptMat(2,2);

figure;
IMG = imread('000_000595-b_8409599_depth.pgm');
subplot(121); imagesc(IMG); title("Depth image");

row = size(IMG,1);
col = size(IMG,2);
cloud = zeros(row*col,3);
k = 0;

for i=1:col
    for j=1:row
        x = (i-cenX)*double(IMG(j,i))/focX;
        y = (j-cenY)*double(IMG(j,i))/focY;
        if(IMG(j,i)~=0)
            k = k+1;
            cloud(k,:) = [x y IMG(j,i)];
        end
    end
end

subplot(122);
plot3(cloud(:,1),cloud(:,2),cloud(:,3),'r.');
grid on;
title("Point cloud");
disp("Perspective projection matrix:");
disp(DeptMat);
