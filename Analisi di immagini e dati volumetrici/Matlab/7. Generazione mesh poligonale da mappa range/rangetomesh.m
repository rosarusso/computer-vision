function [Triangle, Vertex] = rangetomesh(TRUE,XX,YY,ZZ,dim_I,dim_J, K_X84)
% Script per generare una mesh da una mappa range
% Nota che TRUE prende valori 0 o 1
% applica la X84 per eliminare gli edge lunghi
%
% Umberto Castellani
%
%


% crea una struttura per l'indicizzazione dei vertici su un vettore (che
% non posso ricavare dalla mappa range perch� ci sono i buchi)
v_index=ones(dim_I,dim_J);

Vertex=zeros(dim_I*dim_J,3);
ind=0;
for i=1:dim_I
    for j=1:dim_J
        if(TRUE(i,j)==1)
            ind=ind+1;
            v_index(i,j)=ind;
            Vertex(ind,:)=[XX(i,j), YY(i,j), ZZ(i,j)];
        end
    end
end
Vertex=Vertex(1:ind,:);
% Criterio per calcolare il valore di soglia sulle distanze automatico:
res=zeros((dim_I*dim_J)*8,1);
ind=0;
for i=3:dim_I
    for j=3:dim_J
        %controllo gli 8 vicini a nord-est
        if(TRUE(i,j)==1)
            v1=[XX(i,j), YY(i,j), ZZ(i,j)];
            for h=0:2
            for k=0:2
                if(TRUE(i-h,j-k)==1)
                    ind=ind+1;
                    v2=[XX(i-h,j-k), YY(i-h,j-k), ZZ(i-h,j-k)];
                    res(ind)=norm(v1-v2);
                end
            end
            end
        end % if true(i,j)
    end
end
res=res(1:ind);
% Applico la regola X84 che dice:
% |x_i-MED|<k*MAD
MED_x84 = median(res);
MAD_x84 = K_X84 * median(abs(res-MED_x84));
% Poich� sono sicuro che tutti i valori di x_i sono positivi
% questo valore verr� confrontato con tutti i residui dei vertici adiacenti
TX_X84=MED_x84+MAD_x84;


Triangle=zeros((dim_I*dim_J)*8,3);
ind=0;
% creo i triangoli (escludo i bordi, che vengono coinvolti dopo
hw                      =   waitbar(0,'mesh generation ...');
for i=3:dim_I-3
    waitbar(i/(dim_I-3),hw)
    for j=3:dim_J-3
        % controllo che il punto esista:
        if (TRUE(i,j)==1)
            % Devo generare i quattro vertici che determinano 2 triangoli
            ind_central=v_index(i,j);
            central=Vertex(ind_central,:);
            % per ogni vicino devo fare il controllo che sia entro la
            % distaza calcolata con la X84
            %Devo cercare 3 vicini in alto a sinistra:
            %
            %
            % vicino 1: est
            %
            % mi ricavo i punti e le distanze dal pixel centrale
            %
            est1=Vertex(v_index(i,j-1),:);
            diff1=norm(central-est1);
            est2=Vertex(v_index(i,j-2),:);
            diff2=norm(central-est2);
            %
            if (TRUE(i,j-1)==1 && diff1<=TX_X84)
                ind_est=v_index(i,j-1);
                bool_est=1;
                esist_est=1;
            elseif(TRUE(i,j-2)==1 && diff2<=TX_X84)
                ind_est=v_index(i,j-2);
                bool_est=0;
                esist_est=1;
            else
                esist_est=0;
                bool_est=1;
            end
            %
            % vicino 2: nordest
            %
            nordest1=Vertex(v_index(i-1,j-1),:);
            diff1=norm(central-nordest1);
            nordest2=Vertex(v_index(i-1,j-2),:);
            diff2=norm(central-nordest2);
            nordest3=Vertex(v_index(i-2,j-2),:);
            diff3=norm(central-nordest3);
            if (TRUE(i-1,j-1)==1 && diff1<=TX_X84)
                ind_nordest=v_index(i-1,j-1);
                bool_nordest=1;
                esist_nordest=1;
            elseif(bool_est==0 && TRUE(i-1,j-2)==1 && diff2<=TX_X84)
                ind_nordest=v_index(i-1,j-2);
                bool_nordest=1;
                esist_nordest=1; 
            elseif(bool_est==0 && TRUE(i-2,j-2)==1 && diff3<=TX_X84)
                ind_nordest=v_index(i-2,j-2);
                bool_nordest=0;
                esist_nordest=1;
            elseif(bool_est==1 && TRUE(i-2,j-2)==1 && diff3<=TX_X84)    
                ind_nordest=v_index(i-2,j-2);
                bool_nordest=0;
                esist_nordest=1;
            elseif(bool_est==1 && TRUE(i-1,j-2)==1 && diff2<=TX_X84)
                ind_nordest=v_index(i-1,j-2);
                bool_nordest=0;
                esist_nordest=1;
            else
                esist_nordest=0;
                bool_nordest=1;
            end
            %
            % vicino 3: nord
            %
            nord1=Vertex(v_index(i-1,j),:);
            diff1=norm(central-nord1);
            nord2=Vertex(v_index(i-2,j-1),:);
            diff2=norm(central-nord2);
            nord3=Vertex(v_index(i,j-2),:);
            diff3=norm(central-nord3);
            if (TRUE(i-1,j)==1 && diff1<=TX_X84)
                ind_nord=v_index(i-1,j);
                esist_nord=1;
            elseif(TRUE(i-2,j-1)==1 && diff2<=TX_X84)
                ind_nord=v_index(i-2,j-1);
                esist_nord=1; 
            elseif(TRUE(i,j-2)==1 && diff3<=TX_X84)
                ind_nord=v_index(i,j-2);
                esist_nord=1;
            else
                esist_nord=0;
            end
            
        % qui devo creare i 2 triangoli
        % Primo triangolo:
        if (esist_est==1 && esist_nordest==1)
            if (norm(Vertex(ind_central,:)-Vertex(ind_est, :))<=TX_X84 && norm(Vertex(ind_central,:)-Vertex(ind_nordest,:))<=TX_X84 && norm(Vertex(ind_est,:)-Vertex(ind_nordest,:)) <=TX_X84)
            ind=ind+1;
            Triangle(ind,:)=[ind_central ind_est ind_nordest];
            end
        end  
         % Secondo triangolo:
        if (esist_nordest==1 && esist_nord==1)           
            if (norm(Vertex(ind_central,:)-Vertex(ind_nordest, :))<=TX_X84 && norm(Vertex(ind_central,:)-Vertex(ind_nord,:))<=TX_X84 && norm(Vertex(ind_nord,:)-Vertex(ind_nordest,:)) <=TX_X84)
            ind=ind+1;
            Triangle(ind,:)=[ind_central ind_nord ind_nordest];
            end
        end
        % Se il nordest non esiste provo a generarne uno con l'est e il
        % nord
        if (esist_nordest==0 && esist_nord==1 && esist_est==1)
            if (norm(Vertex(ind_central,:)-Vertex(ind_est, :))<=TX_X84 && norm(Vertex(ind_central,:)-Vertex(ind_nord,:))<=TX_X84 && norm(Vertex(ind_nord,:)-Vertex(ind_est,:)) <=TX_X84)
            ind=ind+1;
            Triangle(ind,:)=[ind_central ind_est ind_nord];
            end
        end
        
        end
    end
end
close(hw)
Triangle=Triangle(1:ind,:);
