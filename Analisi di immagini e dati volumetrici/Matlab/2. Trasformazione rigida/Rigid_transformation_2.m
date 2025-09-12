% Trasformazione rigida
% Rosa Russo VR445639

close all
clear all

load Corr3D.mat;

% Cloud points permutation
i = size(data_i,1);
perm = randperm(i);
data = perm(1:ceil(i/3));
data = data_i(data,:);

j = size(model_i,1);
perm = randperm(j);
model = perm(1:ceil(j/3));
model = model_i(model,:);

figure;
plot3(data(:,1),data(:,2),data(:,3),'.b');
hold on;
plot3(model(:,1),model(:,2),model(:,3),'.r');
grid on;

% ICP application
iter = 100;
th = 0;
G = [1 0 0 0
     0 1 0 0
     0 0 1 0
     0 0 0 1];
i = 0;
while(i<iter)
    i = i+1;
    
    % i-th registration
    data = [data_i'; ones(size(data_i,1),1)'];
    datareg = G*data;
    datareg_p = datareg(1:3,:);
    model_p = model';
    
    closest = zeros(size(model));
    min = inf*ones(size(data,1),1);
    % Closest points computation
    for i=1:size(model,1)
        for j=size(datareg_p,1)
            point = norm(model_p(:,i)-datareg_p(:,j));
            
            if(point<min)
                mindist(i) = point;
                closest(i,:) = point;
            end
        end
    end
end