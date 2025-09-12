% test code
close all
clear
clc
DIZ = imread('000_000595-b_8409599_depth.pgm'); %Depth Image z coordinate
figure(1)
imagesc(DIZ);
fx = 575.8;
cx = 319.5;

fy = 575.8;
cy = 239.5;


nrow=size(DIZ, 1); 
ncol=size(DIZ, 2);
cloud=zeros(nrow*ncol, 3);
X=zeros(nrow, ncol);
Y=zeros(nrow, ncol);
TRUE=zeros(nrow, ncol);
k=0;
DephTH=3000;
for i = 1:1:ncol%column
    for j = 1:1:nrow%row9
        X(j,i) = (i - cx) * double(DIZ(j, i)) / fx;
        Y(j,i) = (j - cy) * double(DIZ(j, i)) / fy;
        if (DIZ(j,i)~=0 && DIZ(j,i)<DephTH)
            TRUE(j,i)=1;
            k=k+1;
            cloud(k, :) = [X(j,i) Y(j,i) DIZ(j,i)];
        end
    end
    disp(i);
end

figure(2)
plot3(cloud(:, 1), cloud(:, 2), cloud(:, 3), 'r.');

K_X84=5.2;
[Triangle, Vertex] = rangetomesh(TRUE,X,Y, double(DIZ), nrow, ncol, K_X84);

%Plot del risultato
figure(3)
trisurf(Triangle,Vertex(:,1), Vertex(:,2), Vertex(:,3));
hold on
axis equal

nvertex=size(Vertex,1);
exportMeshToPly(Vertex, Triangle, ones(nvertex,3), 'out2.ply');

