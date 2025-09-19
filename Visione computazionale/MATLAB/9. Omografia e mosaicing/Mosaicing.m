% Homography and mosaicing
% Rosa Russo VR445639

% A homography is a perspective transformation of a plane
% It's a reprojection of a plane from one camera into a different camera
% view, subject to change in the translation (position) and rotation
% (orientation) of the camera.
clear all;
close all;
clc;

%img1=imread('river1.jpg');
%img2=imread('river2.jpg');

%img1 = imread('roofs1.jpg');
%img2 = imread('roofs2.jpg');

%img1 = imread('panorama_image1_big.jpg');
%img2 = imread('panorama_image2_big.jpg');

img1 = imread('city1.jpg');
img2 = imread('city2.jpg');

img1 = rgb2gray(img1);
img2 = rgb2gray(img2);

% Arrays that will contain the chosen points
p1 = [];
p2 = [];

figure(1)
imshow(img1);
figure(2)
imshow(img2);

% Point pick
% We need to choose at least 4 points to calculate a homography between
% the two selected images
for i = 1 : 4 % Number of points to pick
    figure(1);
    hold on
    [x1, y1] = ginput(1); % Pick a point on the first image
    plot(x1, y1, 'r*');
    p1 = [p1; x1, y1];

    figure(2);
    hold on
    [x2, y2] = ginput(1); % Pick a point on the second image
    plot(x2, y2, 'g*'); 
    p2 = [p2; x2, y2];
end

% Complex conjugate transpose.
% It is equal to the transpose if the matrix has non complex numbers.
p1 = p1';
p2 = p2';

% Show selected points
figure;
subplot(121);
imshow(img1);
hold on;
scatter(p1(1,:), p1(2,:), '*r');

subplot(122);
imshow(img2);
hold on;
scatter(p2(1,:), p2(2,:), '*g');

% Homography matrix
% Two images of the same planar surface in space are related by a homography
A = [];
for i = 1 : size(p1,2)
    pi_1 = [p1(:,i);1]; % Homogeneous coordinates -> First image
    pi_1_T = pi_1';
    pi_2 = [p2(:,i);1]; % Homogeneous coordinates -> Second image
    % Outer product
    outer_product = [0 -pi_2(3) pi_2(2); pi_2(3) 0 -pi_2(1); ...
                    -pi_2(2) pi_2(1) 0];
    % Kronecker product
    kronecker = kron(pi_1_T, outer_product);
    % Map Kronecker product over the list of the points and the outer product
    % and build the matrix A
    % A will contain the first two linearly independent rows
    A = [A; kronecker(1,:); kronecker(2,:)];
end
% Transformation construction
h = ns(A);
% H calculation ( H will wrap the second image)
H = reshape(h,3,3);
H = H./H(3,3);

% Building the image from H
[img2, bb, alpha] = imwarp(img2,inv(H), 'linear', 'valid');
ind = find(isnan(img2));
img2(ind)=0;

% Images dimension
[h1, w1] = size(img1);
[h2, w2] = size(img2);
% Computing new corners
x_min = min(bb(1), 1);
x_max = max(bb(3), w1);
y_min = min(bb(2), 1);
y_max = max(bb(4), h1);
% New image dimension
h = y_max - y_min + 1;
w = x_max - x_min + 1;
% img1 and img2 tranformed translations w.r.t the new image
x1_translated = 0;
y1_translated = 0;
x2_translated = 0;
y2_translated = 0;
if (bb(1) < 0)
    x1_translated = -bb(1);
else
    x2_translated = bb(1);
end
if (bb(2) < 0)
    y1_translated = -bb(2);
else
    y2_translated = bb(2);
end

% Composing the image
mosaic = zeros(h, w);
mosaic((y1_translated+1):(y1_translated+h1), (x1_translated+1):(x1_translated+w1)) = img1*0.5;
mosaic((y2_translated+1):(y2_translated+h2), (x2_translated+1):(x2_translated+w2)) = mosaic((y2_translated+1):(y2_translated+h2), (x2_translated+1):(x2_translated+w2)) + img2*0.5;

% Displaying composition
figure;
imshow(uint8(mosaic))