% Trasformazione rigida
% Rosa Russo VR445639

clc
clear all
close all

load Corr3D.mat;

% Centroids computation
centroideX=sum(mean(model_i),1);
centroideY=sum(mean(data_i),1);

% Centralized coordinates
Xi=(model_i-centroideX)';
Yi=(data_i-centroideY)';

% Singular value coordinates decomposition
[U,S,V]=svd(Yi*Xi');

% det=1 -> Rotation matrix
I=eye(3);
I(3,3)=det(V*U');
% Rotation and translation matrices
R=V*I*U';
t=centroideX'-R*centroideY';

% Transformation matrix
G=[  R ,  t];

% Data adaptation in order to have a compatiple multiplication between the
% Transformation matrix and the data
data_imod=[data_i'; ones(size(data_i,1),1)'];

% Registered data
data_out=(G*data_imod)';

figure;
subplot(121);
plot3(model_i(:,1), model_i(:,2), model_i(:,3), 'b.');
hold on;
plot3(data_i(:,1), data_i(:,2), data_i(:,3), 'r.');
grid on;
title('Point clouds before registration');

subplot(122);
plot3(model_i(:,1), model_i(:,2), model_i(:,3), 'b.');
hold on;
plot3(data_out(:,1), data_out(:,2), data_out(:,3), 'r.');
grid on;
title('Point clouds after registration');
axis equal;

disp('Rototranslation matrix:');
disp(G);
