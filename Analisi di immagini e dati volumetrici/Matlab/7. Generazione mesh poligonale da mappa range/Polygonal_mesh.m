% Generazione di una nuvola di punti partendo da una range map 
% Rosa Russo VR445639

clc
close all
clear all

% Intrinsic matrix
DeptMat = [575.8, 0.0, 319.5;
           0.0, 575.8, 239.5;
           0.0, 0.0, 1.0];

% Coordinates
cenX = DeptMat(1,3);
cenY = DeptMat(2,3);
focX = DeptMat(1,1);
focY = DeptMat(2,2);

figure;
IMG=imread('000_000595-b_8409599_depth.pgm');
subplot(121);imagesc(IMG);title("depth image");

row=size(IMG,1);
col=size(IMG,2);
cloud=zeros(row*col,3);
x=zeros(row,col);
y=zeros(row,col);
k=0;
TRUE=zeros(row,col);
DephTH=3000;

for i=1:col
    for j=1:row
        x(j,i)=(i-cenX)*double(IMG(j,i))/focX;
        y(j,i)=(j-cenY)*double(IMG(j,i))/focY;
        
        if(IMG(j,i)~=0 && IMG(j,i)<DephTH)
            TRUE(j,i)=1;
            k=k+1;
            cloud(k,:)=[x(j,i) y(j,i) IMG(j,i)];
        end
    end
end

%% From point cloud to mesh
% Points oordinates x,y,z=IMG(i,j)

% X-84
% |x_i-MED|<k*MAD


% Vertex structure
vidx = ones(row,col); % Vertex index
Vertex = zeros(row*col,3); % Vertex coordinates
idx = 0;
IMG = double(IMG);

for i=1:row
    for j=1:col
        if(TRUE(i,j)==1) % Handling the point cloud "holes"
            idx = idx+1;
            vidx(i,j) = idx;
            Vertex(idx,:) = [x(i,j),y(i,j),IMG(i,j)];
        end
    end
end
Vertex = Vertex(1:idx,:);

% Computing the distance between each vertex and the neighborhood
% (North-East quadrant)
% X-84 rule is applied in order to remove long edges

ris = zeros((row*col)*8,1);
idx = 0;
for i=3:row
    for j=3:col
        if(TRUE(i,j)==1)
            v1 = [x(i,j),y(i,j),IMG(i,j)];
            for h=0:2
                for k=0:2
                    if(TRUE(i-h,j-k)==1)
                        idx = idx+1;
                        v2 = [x(i-h,j-k),y(i-h,j-k),IMG(i-h,j-k)];
                        ris(idx) = norm(v1-v2);
                    end
                end
            end
        end
    end
end
ris = ris(1:idx);

MAD = median(abs(ris-median(ris)));
X84 = (5.2*MAD) + median(ris);

Triangle = zeros((row*col)*8,3);
idx = 0;
% Triangles
for i=3:row-3
    for j=3:col-3
        if (TRUE(i,j)==1) % Checking the presence of the point
            % Generate the 4 vertex determining 2 triangles
            ind_central = vidx(i,j);
            central = Vertex(ind_central,:);
            % For each neighbor, check the distance of the vertex (included
            % in the X-84 calculated distance)
            % Looking for three neighbors inside the north-east quadrant
            % 1st neighbor, East
            est1 = Vertex(vidx(i,j-1),:);
            diff1 = norm(central-est1);
            est2 = Vertex(vidx(i,j-2),:);
            diff2 = norm(central-est2);
            if (TRUE(i,j-1)==1 && diff1<=X84)
                ind_est = vidx(i,j-1);
                bool_est = 1;
                esist_est = 1;
            elseif(TRUE(i,j-2)==1 && diff2<=X84)
                ind_est = vidx(i,j-2);
                bool_est = 0;
                esist_est = 1;
            else
                esist_est = 0;
                bool_est = 1;
            end
            % 2nd neighbor, North-East
            nordest1 = Vertex(vidx(i-1,j-1),:);
            diff1 = norm(central-nordest1);
            nordest2 = Vertex(vidx(i-1,j-2),:);
            diff2 = norm(central-nordest2);
            nordest3 = Vertex(vidx(i-2,j-2),:);
            diff3 = norm(central-nordest3);
            if (TRUE(i-1,j-1)==1 && diff1<=X84)
                ind_nordest = vidx(i-1,j-1);
                bool_nordest = 1;
                esist_nordest = 1;
            elseif(bool_est==0 && TRUE(i-1,j-2)==1 && diff2<=X84)
                ind_nordest = vidx(i-1,j-2);
                bool_nordest = 1;
                esist_nordest = 1; 
            elseif(bool_est==0 && TRUE(i-2,j-2)==1 && diff3<=X84)
                ind_nordest = vidx(i-2,j-2);
                bool_nordest = 0;
                esist_nordest = 1;
            elseif(bool_est==1 && TRUE(i-2,j-2)==1 && diff3<=X84)    
                ind_nordest = vidx(i-2,j-2);
                bool_nordest = 0;
                esist_nordest = 1;
            elseif(bool_est==1 && TRUE(i-1,j-2)==1 && diff2<=X84)
                ind_nordest = vidx(i-1,j-2);
                bool_nordest = 0;
                esist_nordest = 1;
            else
                esist_nordest = 0;
                bool_nordest = 1;
            end
            % 3rd neighbor, North
            nord1 = Vertex(vidx(i-1,j),:);
            diff1 = norm(central-nord1);
            nord2 = Vertex(vidx(i-2,j-1),:);
            diff2 = norm(central-nord2);
            nord3 = Vertex(vidx(i,j-2),:);
            diff3 = norm(central-nord3);
            if (TRUE(i-1,j)==1 && diff1<=X84)
                ind_nord = vidx(i-1,j);
                esist_nord = 1;
            elseif(TRUE(i-2,j-1)==1 && diff2<=X84)
                ind_nord = vidx(i-2,j-1);
                esist_nord = 1; 
            elseif(TRUE(i,j-2)==1 && diff3<=X84)
                ind_nord = vidx(i,j-2);
                esist_nord = 1;
            else
                esist_nord = 0;
            end
            
            % First triangle
            if (esist_est==1 && esist_nordest==1)
                if (norm(Vertex(ind_central,:)-Vertex(ind_est, :))<=X84 && norm(Vertex(ind_central,:)-Vertex(ind_nordest,:))<=X84 && norm(Vertex(ind_est,:)-Vertex(ind_nordest,:)) <=X84)
                idx = idx+1;
                Triangle(idx,:) = [ind_central ind_est ind_nordest];
                end
            end  
             % Second triangle
            if (esist_nordest==1 && esist_nord==1)           
                if (norm(Vertex(ind_central,:)-Vertex(ind_nordest, :))<=X84 && norm(Vertex(ind_central,:)-Vertex(ind_nord,:))<=X84 && norm(Vertex(ind_nord,:)-Vertex(ind_nordest,:)) <=X84)
                idx = idx+1;
                Triangle(idx,:) = [ind_central ind_nord ind_nordest];
                end
            end
            % In case Nosth-East doesn't exists, generate a triangle with
            % East and North quadrants
            if (esist_nordest==0 && esist_nord==1 && esist_est==1)
                if (norm(Vertex(ind_central,:)-Vertex(ind_est, :))<=X84 && norm(Vertex(ind_central,:)-Vertex(ind_nord,:))<=X84 && norm(Vertex(ind_nord,:)-Vertex(ind_est,:)) <=X84)
                idx = idx+1;
                Triangle(idx,:) = [ind_central ind_est ind_nord];
                end
            end
        end
    end
end

Triangle = Triangle(1:idx,:);

figure(3)
trisurf(Triangle,Vertex(:,1), Vertex(:,2), Vertex(:,3));
hold on
axis equal