% Calibration: direct method
% Rosa Russo VR445639

clear all
close all
clc
% Picture of the calibration object
img1 = imread('Image1.jpg');
% Calibration object coordinates (z, x, y)
Mi = [0 0 0 1
      0 0 8.3 1
      5.3 0 8.3 1
      5.3 0 0 1
      0 14.8 8.3 1
      5.3 14.8 8.3 1];
% '1' -> cartesian coordinates
% Show the image
fig1 = figure(1);
imshow(img1);
hold on;

% Collect points: manually select the object vertices
mi = ones(length(Mi), 3);
for i = 1:length(Mi)
    [clickX,clickY] = ginput(1);
    mi(i, :) = [clickX, clickY, 1];
    scatter(clickX, clickY, 'g', '+'); % Plots the coordinates, using '+'
    text(clickX, clickY, strcat('. ',num2str(i)), 'Color', 'g'); % strcat concatenates strings horizontally
end

% Calibration:
A = zeros(2*length(mi), 12);
for i = 1:length(Mi)
    ax = [0, -mi(i,3), mi(i,2);
          mi(i,3), 0, -mi(i,1);
         -mi(i,2), mi(i,1), 0];      
     % Note that Mi is M'
     KRO = kron(ax, Mi(i,:));
     A(2*i-1:2*i, :) = KRO(1:2, :);
end
 [~,~,V] = svd(A,'econ');
 vecP = V(:,end);

P = reshape(vecP, 4, 3).'; % P represents the Perspective Projection Matrix

m_reproj = zeros(length(Mi), 3);
for i = 1:length(Mi)   
    mcurrent = P*(Mi(i,:).');
    m_reproj(i, :) = mcurrent.'/mcurrent(3,1);
end
hold on
for i = 1:length(Mi)
    scatter(m_reproj(i,1), m_reproj(i,2), 'r');
    scatter(mi(i,1), mi(i,2), 'g', '+');
    text((m_reproj(i,1)), (m_reproj(i,2)), strcat('.',num2str(i)));
end
img_i = 1;
% Change the name of PPM here
save(strcat('Calib_direct_', num2str(img_i),'.mat'), 'P');
