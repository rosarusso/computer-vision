% Epipolar lines
% Rosa Russo VR445639

clc
clear all
close all

PPM(:,:,1)=load('Data/angel1.ppm');
PPM(:,:,2)=load('Data/angel2.ppm');
% PPM(:,:,1)=load('left01.ppm');
% PPM(:,:,2)=load('right01.ppm');

% Fundamental matrix F calculation
[F,ep_left,ep_right] = fundamental(PPM(1:3,:,1),PPM(1:3,:,2));

imgleft = imread('Data/angel1.JPG');
imgright = imread('Data/angel2.JPG');
% imgleft = imread('left01.jpg');
% imgright = imread('right01.jpg');

[~,n,~] = size(imgleft);

figure(1);
imshow(imgleft);
i = 3;
title(strcat('Click ', num2str(i),' points'));
axis image;

% Resulting image, containing the epipolar lines corresponding to the
% selected points
figure(2);
imshow(imgright);
title('Corresponding Epipolar Lines');
axis image;
colors = {'yellow', 'green', 'red'};

% Point pick and epipolar line calculation
for i = 1 : i % points to select
    figure(1);
    [x_left, y_left] = ginput(1);
    hold on;
    plot(x_left, y_left,'rx');

    % Finding the epipolar line on the right image:
    p_left = [x_left; y_left; 1];
    % Points needed in order to estimate the line parameters
    p_right = F*p_left;

    epipolar_x = [1;n];
    % Epipolars y from x and point projection on the right image
    epipolar_y = (-p_right(3)-p_right(1)*epipolar_x)/p_right(2);

    figure(2);
    hold on;
    line(epipolar_x, epipolar_y, 'Color', colors{mod(i-1, 3)+1});
end
