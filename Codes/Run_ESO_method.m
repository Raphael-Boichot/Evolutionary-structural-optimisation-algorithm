clear;
clc;
close all;
rng('shuffle', 'twister')
format long
mkdir('Figure');
mkdir('Topology');

%****Thermal properties****************************************************
high_conductivity = 10;         %conductivity of the draining material
low_conductivity = 1;           %conductivity of the heating matter
heat_sink_temperature = 298;    %self explanatory
x_step = 0.001;                 %size of x/x cells
p_vol=1e6;                      %surface of volume power
filling_ratio=0.3;              %ratio of conductive matter on the surface
starting_image='50x100.bmp';    %self explanatory
max_rank=10;                    %maximum rank for exchange
%**************************************************************************

disp('Trying to restart from previous run if any...')
folder=dir('Topology/*.png');
last_valid_file=length(folder);

if not(last_valid_file==0)
    flag=0;
    disp('Previous run detected, loading last known topology...')
    restart_topology=imread(['Topology/',folder(last_valid_file).name]);
    [~,width,~]=size(restart_topology);
    initial_boundary_conditions=[restart_topology(:,1:width/2,:),restart_topology(:,width,:)];
else
    flag=1;
    disp('No preceding run detected, starting from a random topology...')
    initial_boundary_conditions=imread(starting_image);
end

[height,width,profondeur]=size(initial_boundary_conditions);
number_of_images = max([height,width]);
boundary_conditions = zeros(height,width);
history_map=zeros(height,width);
history_map(1,1)=1;
%conversion of an image to boundary conditions
non_conductive_cells=0;
conductive_cells=0;
for k = 1:1:height
    for l = 1:1:width
        red = initial_boundary_conditions(k,l,1);
        green = initial_boundary_conditions(k,l,2);
        blue = initial_boundary_conditions(k,l,3);
        if (red == 255) && (green == 255) && (blue == 255)
            pixel = low_conductivity;
            non_conductive_cells=non_conductive_cells+1;
        end
        if (red == 127) && (green == 127) && (blue == 127)
            pixel = -2;
        end
        if (red == 0) && (green == 0) && (blue == 255)
            pixel = -3;
        end
        if (red == 0) && (green == 0) && (blue == 0)
            pixel = high_conductivity;
            conductive_cells=conductive_cells+1;
        end
        boundary_conditions(k, l) = pixel;
    end
end

if flag==1
    disp('Filling blank image with conductive pixels...')
    number_conductive_cells=ceil(non_conductive_cells*filling_ratio);
    boundary_conditions=init_image(boundary_conditions,number_conductive_cells, low_conductivity, high_conductivity);
end

disp('Starting the ESO algorithm...');

%variable pre-allocation
temp=ones(height,width).*heat_sink_temperature;
boundary_output=initial_boundary_conditions;
affichage=zeros(1,4);
m=last_valid_file;
u=0;
figure('Position',[100 100 600 600]);
local_rank=max_rank;

while max(max(history_map))<50
    tic
    m=m+1;
    disp(' ');
    disp(['---------Epoch: ',num2str(m),'---------']);
    disp('Applying ESO algorithm...');
    [boundary_conditions,growth,etching] = fun_ESO_algorithm(boundary_conditions,high_conductivity,low_conductivity,heat_sink_temperature,x_step,p_vol, max_rank, local_rank);
    [distance,somme_entropie, entropie, border_variance,variance, moyenne_temp,t_max,temp,grad,variance_grad]=finite_temp_direct_sparse(high_conductivity,low_conductivity,heat_sink_temperature,x_step,p_vol,boundary_conditions);
    history_tmax(m-last_valid_file)=t_max;
    
    for k = 1:1:height
        for l = 1:1:width
            pixel = boundary_conditions(k, l) ;
            if pixel == low_conductivity
                red = 255;
                green = 255;
                blue = 255;
            end
            if pixel == -2
                red = 127;
                green = 127 ;
                blue = 127 ;
            end
            if pixel == -3
                red = 0;
                green = 0 ;
                blue = 255;
            end
            if pixel == high_conductivity
                red = 0 ;
                green = 0 ;
                blue = 0;
            end
            boundary_output(k,l,1)=red;
            boundary_output(k,l,2)=green;
            boundary_output(k,l,3)=blue;
        end
    end
    
    boundary_output=uint8(boundary_output);
    mirror=fliplr(boundary_output(1:height,1:width-1,:));
    mirror2=fliplr(mirror);
    arbre=[mirror2,mirror];
    figure(1)
    subplot(2,4,1:2);
    
    if (m-last_valid_file)>2
        for i=2:1:length(history_tmax)
            residuals(i-1)=(history_tmax(i)-history_tmax(i-1))/(history_tmax(2)-history_tmax(1));
        end
        if (m-last_valid_file)<100; plot(log10(abs(residuals)),'.r'); end
        if (m-last_valid_file)>=100; plot(log10(abs(residuals(end-98:end))),'.r'); end
    end
    title('log10 residuals');
    
    subplot(2,4,3:4);
    imagesc([mirror2,mirror]);
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
    title('Gradients');
    
    old_max_history=max(max(history_map));
    for i=1:1:local_rank
        history_map(growth(i,1),growth(i,2))=history_map(growth(i,1),growth(i,2))+1;
        history_map(etching(i,1),etching(i,2))=history_map(etching(i,1),etching(i,2))+1;
    end
    if max(max(history_map))>old_max_history
        local_rank=local_rank-1;
        if local_rank<1
            local_rank=1;
        end
    end
    
    subplot(2,4,8);
    imagesc(log10(history_map));
    title('History map');
    
    disp(['Maximal temperature: ',num2str(history_tmax(m-last_valid_file))]);
    initial_boundary_conditions=boundary_output;
    colormap jet
    pause(0.01);
    
    saveas(gcf,['Figure_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'.png']);
    imwrite(arbre,['Topology_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'.png']);
    saveas(gcf,['Figure/Figure_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'_',num2str(m,'%06.f'),'.png']);
    imwrite(arbre,['Topology/Topology_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'_',num2str(m,'%06.f'),'.png']);
    disp(['Max redunding moves: ',num2str(max(max(history_map)))]);
    disp(['Cells allowed for swapping: ',num2str(local_rank)])
    toc
end
disp('Converged !');
