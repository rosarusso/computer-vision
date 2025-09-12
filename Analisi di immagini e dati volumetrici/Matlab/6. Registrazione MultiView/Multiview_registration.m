% Multiview registration
% Rosa Russo VR445639

clc
close all
clear all

G = eye(4);

figure;
hold on;
grid on;
colore = rand(23,3);
% Data loading and showing
for i=1:23
    if(i<10)
        name = ['./Bunny/bunny-range-00' num2str(i) '.ply'];
    else
        name = ['./Bunny/bunny-range-0' num2str(i) '.ply'];
    end
    
    mesh = readply(name);
    
    % Sub sampled data
    Subdata{i} = [mesh.verts(1:40:end,1),mesh.verts(1:40:end,2),mesh.verts(1:40:end,3)];
    % Full data
    Fulldata{i} = [mesh.verts(:,1),mesh.verts(:,2),mesh.verts(:,3)];
    
    plot3(Subdata{i}(:,1),Subdata{i}(:,2),Subdata{i}(:,3),'.','Color',colore(i,:));
end
data = Subdata{1};
Regdata{1} = Fulldata{1};

for i=2:23
    model = data;
    data = Subdata{i};
    Gnew = icp84(model,data);
    Gloop = G*Gnew;
    Fulldata{i} = [Fulldata{i}'; ones(size(Fulldata{i},1),1)'];
    datareg = (Gloop*Fulldata{i})';
    G = Gloop;
    Regdata{i} = datareg;
end

figure;
hold on;
grid on;
for i=1:23
    plot3(Regdata{i}(:,1), Regdata{i}(:,2), Regdata{i}(:,3), '.', 'Color',colore(i,:));
end