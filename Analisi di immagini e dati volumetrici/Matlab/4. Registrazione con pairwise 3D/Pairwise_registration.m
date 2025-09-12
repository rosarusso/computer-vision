% Registrazione pairwise 3D
% Rosa Russo VR445639

clc
close all
clear all

% Model loading and subsampling
model = load('a4000007.cnn');
perm_m = randperm(size(model,1));
perm_m = perm_m(1:round(size(model,1)/3));
model = model(perm_m,:);

% Data loading and subsampling
data = load('a4000001.cnn');
perm_d = randperm(size(data,1));
perm_d = perm_d(1:round(size(data,1)/3));
data = data(perm_d,:);

figure;
subplot(121);
title('Pre ICP');
plot3(model(:,1), model(:,2), model(:,3), '.b');
hold on;
plot3(data(:,1), data(:,2), data(:,3), '.r');
grid on;
axis equal;

%% ICP algorithm
threshold = 0.00000001;
G = eye(4);
ris = inf;
risprev = 0;
i = 0;

while((abs(ris-risprev)>threshold) && i<200)
    i = i + 1;
    risprev = ris;
    G = G(1:3,:);
    % The new transformation matrix is applied to the data
    datamod = [data'; ones(size(data,1),1)'];
    dataReg = (G*datamod)';
    
    closest = zeros(size(data));
    mindist = inf*ones(size(data,1),1);
    % Closest points computation
    for j=1:size(dataReg,1)
        for z=1:size(model,1)
            d = norm(model(z,:) - dataReg(j,:));
            if(d<mindist(j))
                mindist(j) = d;
                closest(j,:) = model(z,:);
            end
        end
    end
    
    ris = mean(mindist);
    
    % Orthogonal procustian problem
    % Centroids
    centroideX = sum(mean(closest),1);
    centroideY = sum(mean(dataReg),1);

    % Centralized coordinates
    Xi = (closest-centroideX)';
    Yi = (dataReg-centroideY)';

    % Singular value centralized coordinates decomposition
    [U,S,V]=svd(Yi*Xi');

    % Matrix that allows to obtain rotation matrix with det=1
    I=eye(3);
    I(3,3)=det(V*U');

    % Rotation and translation matrix
    R=V*I*U';
    t=centroideX'-R*centroideY';

    % Transformation matrix
    Gnew=[  R ,  t;
           0 0 0 1];
    G=[G; 0 0 0 1];
    G=Gnew*G;
end

% applico la matrice di trasformazione finale ottenuta ai dati iniziali
final_data = (G*datamod)';

subplot(122);
title('Post-ICP');
plot3(model(:,1), model(:,2), model(:,3), '.b');
hold on
plot3(final_data(:,1), final_data(:,2), final_data(:,3), '.r');
hold on
grid on
axis equal
