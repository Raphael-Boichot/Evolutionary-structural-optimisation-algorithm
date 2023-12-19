function [condi_limites,growth,etching] = fun_ESO_algorithm(condi_limites,cond_haute,cond_basse,temp_puits,pas_x,p_vol)
[hauteur,largeur,profondeur]=size(condi_limites);
%******************************************Début de l'automate cellulaire**
%identification des pixels à faire pousser
rng('shuffle', 'twister')
grow=zeros(hauteur,largeur);
for k = 1:1:hauteur;
    for l = 1:1:largeur;
        if condi_limites(k,l)==cond_haute;
            if condi_limites(k+1,l)==cond_basse; grow(k+1,l)=1; end;
            if condi_limites(k-1,l)==cond_basse; grow(k-1,l)=1; end;
            if condi_limites(k,l+1)==cond_basse; grow(k,l+1)=1; end;
            if condi_limites(k,l-1)==cond_basse; grow(k,l-1)=1; end;
        end;
    end;
end;

m=0;
for k = 1:1:hauteur;
    for l = 1:1:largeur;
        if grow(k,l)==1;
            m=m+1;
            grow_pos(m,1)=k;
            grow_pos(m,2)=l;
        end;
    end;
end;

%******************************************Début de l'automate cellulaire**
%identification des pixels à faire crever
etch=zeros(hauteur,largeur);
for k = 1:1:hauteur;
    for l = 1:1:largeur;
        if condi_limites(k,l)==cond_basse;
            if condi_limites(k+1,l)==cond_haute; etch(k+1,l)=1; end;
            if condi_limites(k-1,l)==cond_haute; etch(k-1,l)=1; end;
            if condi_limites(k,l+1)==cond_haute; etch(k,l+1)=1; end;
            if condi_limites(k,l-1)==cond_haute; etch(k,l-1)=1; end;
        end;
    end;
end;

m=0;
for k = 1:1:hauteur;
    for l = 1:1:largeur;
        if etch(k,l)==1;
            m=m+1;
            etch_pos(m,1)=k;
            etch_pos(m,2)=l;
        end;
    end;
end;



%identification des pixels à faire pousser
condi_limites_temp=zeros(hauteur,largeur,length(grow_pos));
for m=1:1:length(grow_pos);
    condi_temp=condi_limites;
    condi_temp(grow_pos(m,1),grow_pos(m,2))=cond_haute;
    condi_limites_temp(:,:,m)=condi_temp;
end
parfor m=1:1:length(grow_pos);
    [somme_entropie, entropie, border_variance,variance, moyenne_temp,temp_max,temp,grad, variance_grad]=finite_temp_direct_sparse(cond_haute,cond_basse,temp_puits,pas_x,p_vol,condi_limites_temp(:,:,m));
    grow_pos(m,3)=temp_max;
end

grow_pos=sortrows(grow_pos,3);

%identification des pixels à faire crever
condi_limites_temp=zeros(hauteur,largeur,length(etch_pos));
for m=1:1:length(etch_pos);
    condi_temp=condi_limites;
    condi_temp(etch_pos(m,1),etch_pos(m,2))=cond_basse;
    condi_limites_temp(:,:,m)=condi_temp;
end
parfor m=1:1:length(etch_pos);
    [somme_entropie, entropie, border_variance,variance, moyenne_temp,temp_max,temp,grad, variance_grad]=finite_temp_direct_sparse(cond_haute,cond_basse,temp_puits,pas_x,p_vol,condi_limites_temp(:,:,m));
    etch_pos(m,3)=temp_max;
end
etch_pos=sortrows(etch_pos,3);

growth=grow_pos(ceil(5*rand),:);
etching=etch_pos(ceil(5*rand),:);

%Echange des cellules
condi_limites(etching(1,1),etching(1,2))=cond_basse;
condi_limites(growth(1,1),growth(1,2))=cond_haute;
%******************************************Fin de l'automate cellulaire****
