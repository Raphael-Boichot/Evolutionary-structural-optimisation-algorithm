clear;
clc;
close all;
format long
rng('shuffle', 'twister')
mkdir('Figure');
mkdir('Topology');

%****Thermal properties******************************
high_conductivity = 10;
low_conductivity = 1;
heat_sink_temperature = 300;
x_step = 0.001;
p_vol=1e5;
filling_ratio=0.3;
starting_image='50x100.bmp';

disp('Trying to restart from last good image...')
folder=dir('Topology/*.png');
last_valid_file=length(folder);

if not(last_valid_file==0)
    flag=0;
    restart_topology=imread(['Topology/',folder(last_valid_file).name]);
    [~,largeur,~]=size(restart_topology);
    condi_limites_1=[restart_topology(:,1:largeur/2,:),restart_topology(:,largeur,:)];
else
    flag=1;
    disp('No preceding run detected, starting from scratch...')
    condi_limites_1=imread(starting_image);
end

[hauteur,largeur,profondeur]=size(condi_limites_1);
nombre_images = max([hauteur,largeur]);
condi_limites = zeros(hauteur,largeur);
automate=zeros(hauteur,largeur);
condi_temp = zeros(hauteur,largeur);
f_temp=zeros(nombre_images);
f_iter=1:1:nombre_images;

%conversion of an image to boundary conditions
pixels_blancs=0;
pixels_noirs=0;
for k = 1:1:hauteur
    for l = 1:1:largeur
        rouge = condi_limites_1(k,l,1);
        vert = condi_limites_1(k,l,2);
        bleu = condi_limites_1(k,l,3);
        
        if (rouge == 255) && (vert == 255) && (bleu == 255)
            choix = low_conductivity;
            pixels_blancs=pixels_blancs+1;
        end
        if (rouge == 127) && (vert == 127) && (bleu == 127)
            choix = -2;
        end
        if (rouge == 0) && (vert == 0) && (bleu == 255)
            choix = -3;
        end
        if (rouge == 0) && (vert == 0) && (bleu == 0)
            choix = high_conductivity;
            pixels_noirs=pixels_noirs+1;
        end
        
        condi_limites(k, l) = choix;
        
    end
end

if flag==1
disp('Filling image with conductive pixels...')
nombre_pixels_conducteurs=ceil(pixels_blancs*filling_ratio);
condi_limites=init_image(condi_limites,nombre_pixels_conducteurs, low_conductivity, high_conductivity);
end

disp('Starting ESO algorithm...');

%variable pre-allocation
temp=ones(hauteur,largeur).*heat_sink_temperature;
condu_tab=zeros(hauteur,largeur,4);
new_pos_in=zeros(hauteur,largeur);
new_pos_out=zeros(hauteur,largeur);
new_pos2=zeros(hauteur,largeur);
gradients=zeros(hauteur,largeur);
note=zeros(hauteur,largeur);
condi_limites_2=condi_limites_1;
affichage=zeros(1,4);

m=last_valid_file;
u=0;
figure('Position',[100 100 600 600]);

while max(max(automate))<20
    tic
    m=m+1;
    disp(['---------Epoch: ',num2str(m),'---------']);
    disp('Applying ESO algorithm...');
    [condi_limites,growth,etching] = fun_ESO_algorithm(condi_limites,high_conductivity,low_conductivity,heat_sink_temperature,x_step,p_vol);
    [somme_entropie, entropie, border_variance,variance, moyenne_temp,t_max,temp,grad, variance_grad]=finite_temp_direct_sparse(high_conductivity,low_conductivity,heat_sink_temperature,x_step,p_vol,condi_limites);
    t_max_sortie(m-last_valid_file)=t_max;
    for k = 1:1:hauteur
        for l = 1:1:largeur
            
            choix = condi_limites(k, l) ;
            
            if choix == low_conductivity
                rouge = 255;
                vert = 255;
                bleu = 255;
            end
            
            if choix == -2
                rouge = 127;
                vert = 127 ;
                bleu = 127 ;
            end
            if choix == -3
                rouge = 0;
                vert = 0 ;
                bleu = 255;
            end
            if choix == high_conductivity
                rouge = 0 ;
                vert = 0 ;
                bleu = 0;
            end
            
            condi_limites_2(k,l,1)=rouge;
            condi_limites_2(k,l,2)=vert;
            condi_limites_2(k,l,3)=bleu;
            
        end
    end
    
    
    condi_limites_2=uint8(condi_limites_2);
    miroir=fliplr(condi_limites_2(1:hauteur,1:largeur-1,:));
    miroir2=fliplr(miroir);
    arbre=[miroir2,miroir];
    figure(1)
    subplot(2,4,1:2);
    
    
    if (m-last_valid_file)>2
        for i=2:1:length(t_max_sortie)
            residual(i-1)=(t_max_sortie(i)-t_max_sortie(i-1))/(t_max_sortie(2)-t_max_sortie(1));
        end
        
        if (m-last_valid_file)<100; plot(log10(abs(residual)),'.r'); end
        if (m-last_valid_file)>=100; plot(log10(abs(residual(end-98:end))),'.r'); end
        
    end
    
    title('log10 Residuals');
    
    subplot(2,4,3:4);
    imagesc([miroir2,miroir]);
    title('Topology');
    
    subplot(2,4,5);
    plot(variance,'.m');
    imagesc(log10(entropie(2:end-1,2:end-1)));
    title('Log10 Entropy');
    
    subplot(2,4,6);
    imagesc(temp);
    title('Temperature');
    
    subplot(2,4,7);
    imagesc(grad);
    colormap hot
    title('Gradients');
    
    subplot(2,4,8);
    automate(growth(1,1),growth(1,2))=automate(growth(1,1),growth(1,2))+1;
    automate(etching(1,1),etching(1,2))=automate(etching(1,1),etching(1,2))+1;
    imagesc(automate);
    title('Redundancy');
    
    disp(['Maximal temperature: ',num2str(t_max_sortie(m-last_valid_file))]);
    condi_limites_1=condi_limites_2;
    colormap jet
    pause(0.01);
    
    saveas(gcf,['Figure_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'.png']);
    imwrite(arbre,['Topology_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'.png']);
    saveas(gcf,['Figure/Figure_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'_',num2str(m,'%06.f'),'.png']);
    imwrite(arbre,['Topology/Topology_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'_',num2str(m,'%06.f'),'.png']);
    toc
end
disp('Converged !');
