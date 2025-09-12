%
% Script for rigid pose estimation
%
clear all
close all

load Corr3D.mat



figure(1);
plot3(model_i(:,1), model_i(:,2), model_i(:,3), 'b.');
hold on
plot3(data_i(:,1), data_i(:,2), data_i(:,3), 'r.');
grid on
axis equal
%
%
G_out=absolute(model_i, data_i);
data_out = rigid(G_out,data_i);
%
figure(2);
plot3(model_i(:,1), model_i(:,2), model_i(:,3), 'b.');
hold on
plot3(data_out(:,1), data_out(:,2), data_out(:,3), 'r.');
grid on
axis equal
