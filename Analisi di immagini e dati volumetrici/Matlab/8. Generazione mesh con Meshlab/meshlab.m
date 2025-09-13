% Generazione di mesh con Meshlab
% Rosa Russo VR445639
% Load point cloud
ptCloud = pcread('PCregistered.ply');
points = double(full(real(ptCloud.Location)));

alpha_values = linspace(1, 10, 10); % Alpha range

figure('Name', 'AlphaShape Reconstruction with Rotation', 'Position', [200 200 1200 800]);

for i = 1:length(alpha_values)
    alphaRadius = alpha_values(i);
    shp = alphaShape(points, alphaRadius);

    subplot(2,5,i);
    plot(shp);
    axis equal;
    title(['Alpha = ' num2str(alphaRadius)]);
    view(3); % Set 3D view mode
    camorbit(36*(i-1), 0); % Rotate camera horizontally by increments of 36 deg per plot
    drawnow;
end

