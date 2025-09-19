% Rigid pose estimation
% Rosa Russo VR445639

clear
close all
clc

load Corr3D.mat

figure(1);
plot3(model_i(:,1), model_i(:,2), model_i(:,3), 'b.');
hold on
plot3(data_i(:,1), data_i(:,2), data_i(:,3), 'r.');
grid on
axis equal

% Compute rigid transformation between two 3D point sets "model_i" and "data_i"
G_out = absolute(model_i, data_i);
% Apply 3D rigid transformation G to a point set "data_i"
data_out = rigid(G_out,data_i);

figure(2);
plot3(model_i(:,1), model_i(:,2), model_i(:,3), 'b.');
hold on
plot3(data_out(:,1), data_out(:,2), data_out(:,3), 'r.');
title('Compute the rigid transformation between two 3D point sets');
grid on
axis equal
