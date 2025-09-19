% Triangulation
% Rosa Russo VR445639

clear all
close all
clc

load  Calib_direct_1.mat
P1 = P;

load  Calib_direct_2.mat
P2 = P;

% Projection Matrix
PPM = {P1,P2};

I1 = imread('Image1.jpg');
I2 = imread('Image2.jpg');

% Pick 2 points
figure;
imshow(I1);
p1 = ginput(1);
hold on;
plot(p1(1),p1(2),'g*');

p2 = ginput(1);
hold on;
plot(p2(1),p2(2),'g*');
mleft = {p1, p2};

% Pick the same points of the previous image
figure;
imshow(I2);
p1 = ginput(1);
hold on;
plot(p1(1),p1(2),'g*');

p2 = ginput(1);
hold on;
plot(p2(1),p2(2),'g*');
mright = {p1, p2};
m = {mleft,mright};

% Model extimation
M = intersect_base(PPM,m);

% Calculation of the Eulerian distance between the chosen points
d = sqrt(sum( (M(:,1) - M(:,2)) .^2));
close all

% Results
imshow(I1);
hold on, line(cellfun(@(x) x(1), mleft),cellfun(@(x) x(2), mleft), 'Color', 'green')
hold on, text((mleft{1}(1)+mleft{2}(1))/2+10,(mleft{1}(2)+mleft{2}(2))/2+10,num2str(d),'FontSize',16,'Color','g')