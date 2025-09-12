% Registrazione con algoritmo X-84 ICP
% Rosa Russo VR445639

clc
close all
clear all

% Model loading and subsampling 
model = load('a4000007.cnn');
perm_m = randperm(size(model,1));
perm_m = perm_m(1:round(size(model,1)/3));
modelperm = model(perm_m,:);

% Data loading and subsampling
data = load('a4000001.cnn');
outliers = data(1000:1500,:)+randn(501,3).*50+ones(501,3).*300;
data = [data;outliers];
perm_d = randperm(size(data,1));
perm_d = perm_d(1:round(size(data,1)/3));
dataperm = data(perm_d,:);

figure;
subplot(121);
title('Pre ICP');
plot3(model(:,1), model(:,2), model(:,3), '.b');
hold on;
plot3(data(:,1), data(:,2), data(:,3), '.r');
grid on;
axis equal;

% ICP algorithm
treshold = 0.00000001;
G = eye(4);
ris = inf;
risprev = 0;
i = 0;

while((abs(ris-risprev)>treshold) && i<200)
    i = i+1;
    risprev = ris;
    G = G(1:3,:);
    
    datamod=[dataperm'; ones(size(dataperm,1),1)'];
    dataReg=(G*datamod)';
    
    closest=zeros(size(dataperm));
    mindist = inf*ones(size(dataperm,1),1);
    % Closest points
    for j=1:size(dataReg,1)
        for z=1:size(modelperm,1)
            d=norm(modelperm(z,:)-dataReg(j,:));
            if(d<mindist(j))
                mindist(j)=d;
                closest(j,:)=modelperm(z,:);
            end
        end
    end
    
    % X-84 rule application
    MAD = median(abs(mindist-median(mindist)));
    out = 5.2*MAD;
    
    % Outliers/inliers point checking
    scarto = find(abs(mindist-median(mindist))>out);
    closest(scarto,:) = NaN;
    inliers = find(abs(mindist-median(mindist))<=out);
    ris = mean(mindist(inliers));
    
    % Outliers elimination
    delete = find(~isnan(closest));
    closest = reshape(closest(delete),length(closest(delete))/3,3);
    dataReg = reshape(dataReg(delete),length(dataReg(delete))/3,3);
    
    % Orthogonal Procustian problem
    centroideX = sum(mean(closest),1);
    centroideY = sum(mean(dataReg),1);

    % Centralized coordinates
    Xi = (closest-centroideX)';
    Yi = (dataReg-centroideY)';

    % Singular value centralized coordinates decomposition
    [U,S,V] = svd(Yi*Xi');

    % det=1 -> Rotation matrix
    I = eye(3);
    I(3,3) = det(V*U');

    % Rotation and translation matrix
    R = V*I*U';
    t = centroideX'-R*centroideY';

    % Transformation matrix
    Gnew = [  R ,  t;
             0 0 0 1];
    G = [G; 0 0 0 1];
    G = Gnew*G;
end

% Transformation matrix applied to initial data
final_data = (G*datamod)';

subplot(122);
title('Post-ICP');
plot3(modelperm(:,1), modelperm(:,2), modelperm(:,3), '.b');
hold on
plot3(final_data(:,1), final_data(:,2), final_data(:,3), '.r');
hold on
grid on
axis equal
