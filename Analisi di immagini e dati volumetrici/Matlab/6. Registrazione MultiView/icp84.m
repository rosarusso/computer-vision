function G = icp84(model,data)
    soglia=0.00000001;

    G=eye(4);
    ris=inf;
    risprev=0;
    i=0;

    while((abs(ris-risprev)>soglia) && i<200)
        i=i+1;
        risprev=ris;
        G=G(1:3,:);
        %applico la nuova matrice di trasformazione ai dati
        datamod=[data'; ones(size(data,1),1)'];
        dataReg=(G*datamod)';

        closest=zeros(size(data));
        mindist = inf*ones(size(data,1),1);
        %trovo i closest point
        for j=1:size(dataReg,1)
            for z=1:size(model,1)
                d=norm(model(z,:)-dataReg(j,:));
                if(d<mindist(j))
                    mindist(j)=d;
                    closest(j,:)=model(z,:);
                end
            end
        end

        %applicazione della regola X84
        MAD=median(abs(mindist-median(mindist)));
        out=5.2*MAD;

        %verifica dei punti che sono o meno outliers
        scarto=find(abs(mindist-median(mindist))>out);
        closest(scarto,:)=NaN;

        inliers=find(abs(mindist-median(mindist))<=out);
        ris=mean(mindist(inliers));

        %prima di risolvere il problema procustiano ortogonale
        %scarto gli outliers
        eliminare = find(~isnan(closest));
        closest = reshape(closest(eliminare),length(closest(eliminare))/3,3);
        dataReg = reshape(dataReg(eliminare),length(dataReg(eliminare))/3,3);

        %%risolvo il problema procustiano ortogonale
        %calcolo i centroidi
        centroideX=sum(mean(closest),1);
        centroideY=sum(mean(dataReg),1);

        % coordinate centralizzate
        Xi=(closest-centroideX)';
        Yi=(dataReg-centroideY)';

        %scompongo le coord centralizzate per singular value decomposition
        [U,S,V]=svd(Yi*Xi');

        %creo la matrice per far venire la matrice di rotazione con det=1
        I=eye(3);
        I(3,3)=det(V*U');

        %creo matrice di rotazione e traslazione
        R=V*I*U';
        t=centroideX'-R*centroideY';

        %creo la matrice di trasformazione
        Gnew=[  R ,  t;
               0 0 0 1];
        G=[G; 0 0 0 1];
        G=Gnew*G;
    end
end

