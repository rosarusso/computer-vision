% Simulation of the animation
% Rosa Russo VR445639

clear all
close all
clc
img = imread('Image1.jpg');
load('Calib_direct_1.mat');

% Animation 1
% Projection and animation of a 3D point
imshow(img)
hold on
% Define 3D points as rows with homogeneous coordinates [X Y Z 1]
Mi = [0 0 0 1
      0 0 8.3 1
      5.3 0 8.3 1
      5.3 0 0 1
      0 14.8 8.3 1
      5.3 14.8 8.3 1];
% Starting point: corner of the box
point = [5.3 0 8.3 1]';
% Project 3D point into 2D image coordinates using camera matrix P
m = P*point;
m = m / m(3); % Homogeneous normalization, last coordinate is 1
% Plot the projected point as a red circle
pl = plot(m(1), m(2), 'ro', 'MarkerSize', 5, 'LineWidth', 2);
for i = 2 : 51
    pause(0.2);
    % Move from 0 to 14.8
    point = [5.3, 14.8 * (i - 1) / 50, 8.3, 1]';
    m = P * point;
    m = m / m(3);
    % Update plotted point position on the image
    set(pl, 'XData', m(1), 'YData', m(2))
end

% Animation 2
% Projection of a bouncing sphere
img = imread('Image2.jpg');
load('Calib_direct_2.mat');

figure(2);
[X,Y,Z] = sphere;

% filename = 'animation_sphere.gif'; % Output filename

for i = 1:3
    % Projection of the sphere translated by center(t)
    for t = 0:0.1:3
        imshow(img);
        hold on;
        %[x,y,z] = center(t); % Sphere center
        % Define sphere center position based on time t
        % x moves linearly, y constant, z follows quadratic curve
        x = -0.47*t+2.9;
        y = 9.3;
        z = 2*t^2-1.6400*t-1.8000;
        Sph = [1.8*X(:).'+x; 1.8*Y(:).'+y; 1.8*Z(:).'+z; ones(1,numel(X))];
        % Normalize homogeneous coordinates
        sph = normalize(P*Sph);
        % Plot the projected sphere points as scatter plot
        scatter(sph(1,:),sph(2,:));
        hold off;
        pause(0.01);
        %drawnow;
        % % Capture frame and write to GIF
        % frame = getframe(gcf);
        % im = frame2im(frame);
        % [imind, cm] = rgb2ind(im, 256);
        %
        % if i == 1 && t == 0
        %     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
        % else
        %     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
        % end
        pause(0.01);
    end
    for t = 3:-0.1:0
        imshow(img);
        hold on;
        [x,y,z] = center(t);
        Sph = [1.8*X(:).'+x; 1.8*Y(:).'+y; 1.8*Z(:).'+z; ones(1,numel(X))];
        sph = normalize(P*Sph);
        scatter(sph(1,:),sph(2,:));

        hold off;
        pause(0.005);

        %  % Capture frame and write to GIF
        % frame = getframe(gcf);
        % im = frame2im(frame);
        % [imind, cm] = rgb2ind(im, 256);
        %
        % if i == 1 && t == 0
        %     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
        % else
        %     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
        % end
        pause(0.01);
    end
end



% Sphere center
function [x,y,z] = center(t)
   x = -0.47*t+2.9;
   y = 9.3;
   z = 2*t^2-1.6400*t-1.8000;
end




% Animation 3
img = imread('Image1.jpg');
load('Calib_direct_1.mat');

figure(3);
imshow(img);
hold on;

pos = [4.5, 4, 4];  % Sphere center coordinates [x,y,z]
vel = [0, 0.3, 0]; % Velocity vector, moving in y-axis only
raggio = 1.8;       % Sphere radius

% 3D bounding box limits only in y-axis
boxLimits.y = [0+raggio/2 14.8-raggio/2]; % Sphere center coordinates [x,y,z]

[X, Y, Z] = sphere(30);

for frame = 1:300
    % Update sphere position based on velocity
    pos = pos + vel;

    % Check for collisions with box limits to simulate bounce
    if pos(2) - raggio < boxLimits.y(1)
        vel(2) = abs(vel(2)); % Bounce
    elseif pos(2) + raggio > boxLimits.y(2)
        vel(2) = -abs(vel(2)); % Bounce
    end

    % Generate 3D coordinates of sphere surface translated to current position
    Sph = [raggio*X(:).' + pos(1);
           raggio*Y(:).' + pos(2);
           raggio*Z(:).' + pos(3);
           ones(1,numel(X))];

    proj = P * Sph;
    proj = proj ./ proj(3,:); % Normalize homogeneous coordinates

    imshow(img);
    hold on;
    scatter(proj(1,:), proj(2,:), 10, 'r', 'filled'); % Draw projected sphere as red dots
    hold off;

    pause(0.02);
end

% Animation 4
img = imread('Image1.jpg');
load('Calib_direct_1.mat');

figure(4);
imshow(img);
hold on;

% Define 3D vertices of the box in homogeneous coordinates
Mi = [ ...
    0 0 0 1;        % 1: bottom-back-left
    0 0 8.3 1;      % 2: bottom-front-left
    5.3 0 8.3 1;    % 3: bottom-front-right
    5.3 0 0 1;      % 4: bottom-back-right
    0 14.8 8.3 1;   % 5: top-front-left
    5.3 14.8 8.3 1; % 6: top-front-right
    0 14.8 0 1;     % 7: top-back-left
    5.3 14.8 0 1];  % 8: top-back-right

edges = [ ...
    4 1; 1 2; 2 3; 3 4; ... % Left face edges
    4 8; 8 6; 6 3; ... % Top face edges, and connection from top to bottom
    3 2; 2 5; 5 6];    % Front edges connecting bottom and top faces

colors = ['r' 'g' 'b' 'm' 'c' 'y' 'k'];

pl = [];
framesPerEdge = 50;
% filename = 'animation_edges.gif'; % Output filename

for e = 1:size(edges,1)
    startPt = Mi(edges(e,1), :)';
    endPt = Mi(edges(e,2), :)';
    cIndex = mod(e-1, length(colors)) + 1;
    color = colors(cIndex);

    for k = 0:framesPerEdge
        point = startPt * (1 - k/framesPerEdge) + endPt * (k/framesPerEdge);
        m = P * point;
        m = m / m(3);

        if isempty(pl)
            if length(color) == 3
                pl = plot(m(1), m(2), 'o', 'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', color, 'MarkerFaceColor', color);
            else
                pl = plot(m(1), m(2), [color 'o'], 'MarkerSize', 7, 'LineWidth', 2);
            end
        else
            set(pl, 'XData', m(1), 'YData', m(2), 'MarkerEdgeColor', color, 'MarkerFaceColor', color);
        end

        drawnow;
        % % Capture the plot as an image
        % frame = getframe(gcf);
        % im = frame2im(frame);
        % [imind, cm] = rgb2ind(im, 256);
        %
        % % Write to the GIF file
        % if e == 1 && k == 0
        %     imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
        % else
        %     imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
        % end
    end
end


