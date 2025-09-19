% Rectification
% Rosa Russo VR445639

close all
clc
stereo = load('Calib_Results_stereo.mat');

imgleft = imread('left.jpg');
imgright = imread('right.jpg');

% Setting parameters

% Left P
P_left = stereo.KK_left*[1 0 0 0;
                         0 1 0 0;
                         0 0 1 0]*eye(4);

% Right P
P_right = stereo.KK_right*[1 0 0 0;
                           0 1 0 0;
                           0 0 1 0]*[stereo.R stereo.T;
                                     0 0 0 1];

Q_left = P_left(:,1:3);
Q_right = P_right(:,1:3);

q_left = P_left(:,4);
q_right = P_right(:,4);

% PS (Product Sum) method for the factorization
[Kl, R_left, t_left] = KRt_fact(P_left);
[Kr, R_right, t_right] = KRt_fact(P_right);

% Centers
c_left = - inv(Q_left)*q_left;
c_right = - inv(Q_right)*q_right;

% Centers axis
v0 = (c_right - c_left);

% Computing new axis
k = R_left(3,:);
v1 = cross(k',v0);
v2 = cross(v0,v1);

% Rotation
R = [v0'/norm(v0);
     v1'/norm(v1);
     v2'/norm(v2)];

% New PS, given the rotation R
Pnleft = Kl * [R -R*c_left];
Pnright = Kr * [R -R*c_right];

% Transformation matrixes
T_left = Pnleft(:,1:3) /Q_left;
T_right = Pnright(:,1:3) /Q_right;

%% IMAGE RECTIFICATION
imgleft_rect = imwarp(imgleft,T_left);
imgright_rect = imwarp(imgright,T_right);

% Fundamental calculation
[F,eleft,eright] = fundamental(Pnleft,Pnright);

[~,n,~] = size(imgleft_rect);

figure(1);
imshow(uint8(imgleft_rect));
title('Click 3 points');
% Resulting image, containing the epipolar lines corresponding to the
% selected points
figure(2);
imshow(uint8(imgright_rect));
title('Corresponding Epipolar Lines');
colors = {'yellow', 'green', [1 0.6 0.1]};

% Point pick and epipolar line calculation
for i = 1 : 3 % 3 points selected
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