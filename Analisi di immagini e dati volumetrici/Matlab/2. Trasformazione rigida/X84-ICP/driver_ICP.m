% Driver per test di ICP

% Author: A. Fusiello, 2003
% revised by U. Castellani 2003
clear 
close all
% Set parameters
srate = 3;

data = load('a4000001.cnn');
% Subsample data
s = size(data,1);
i = randperm(s);
i = i(1:round(s/srate));
data = data(i,:);

model = load('a4000007.cnn');
% Subsample model
s = size(model,1);
i = randperm(s);
i = i(1:round(s/srate));
model = model(i,:);

figure(100);
plot3(model(:,1), model(:,2), model(:,3), '.b');
hold on
plot3(data(:,1), data(:,2), data(:,3), '.r');
hold on
grid on
axis equal

Gi = icp(model,data);
datafinal = rigid(Gi,data);

figure(200);
plot3(model(:,1), model(:,2), model(:,3), '.b');
hold on
plot3(datafinal(:,1), datafinal(:,2), datafinal(:,3), '.r');
hold on
grid on
axis equal
