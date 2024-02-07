clear;
clc;
close all;
format long
rng('shuffle', 'twister')
mkdir('Figure');
mkdir('Topology');

%****Thermal properties**************************
high_conductivity = 10;
low_conductivity = 1;
heat_sink_temperature = 300;
x_step = 0.001;
p_vol=1e5;
filling_ratio=0.3;
%****ESO parameters******************************
rank_cells_to_exchange=5;%maximal rank allowed for exchange
number_cells_allowed_to_move=1;%must be less than or equal to rank_cells_to_exchange
starting_image='50x100.bmp';

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

while max(max(history_map))<30
    tic
    m=m+1;
    disp(' ');
    disp(['---------Epoch: ',num2str(m),'---------']);
    disp('Applying ESO algorithm...');
    [boundary_conditions,growth,etching] = fun_ESO_algorithm(boundary_conditions,high_conductivity,low_conductivity,heat_sink_temperature,x_step,p_vol,rank_cells_to_exchange,number_cells_allowed_to_move);
    [somme_entropie, entropie, border_variance,variance, moyenne_temp,t_max,temp,grad, variance_grad]=finite_temp_direct_sparse(high_conductivity,low_conductivity,heat_sink_temperature,x_step,p_vol,boundary_conditions);
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
    colormap hot
    title('Gradients');
    
    subplot(2,4,8);
    for i=1:1:number_cells_allowed_to_move
        history_map(growth(i,1),growth(i,2))=history_map(growth(i,1),growth(i,2))+1;
        history_map(etching(i,1),etching(i,2))=history_map(etching(i,1),etching(i,2))+1;
    end
    imagesc(history_map);
    title('History');
    
    disp(['Maximal temperature: ',num2str(history_tmax(m-last_valid_file))]);
    initial_boundary_conditions=boundary_output;
    colormap jet
    pause(0.01);
    
    saveas(gcf,['Figure_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'.png']);
    imwrite(arbre,['Topology_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'.png']);
    saveas(gcf,['Figure/Figure_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'_',num2str(m,'%06.f'),'.png']);
    imwrite(arbre,['Topology/Topology_kp_ko_',num2str(high_conductivity),'_phi_',num2str(filling_ratio),'_',num2str(m,'%06.f'),'.png']);
    
    if (m-last_valid_file)>2
        if residuals(end)<0
            number_cells_allowed_to_move=number_cells_allowed_to_move-1;
            if number_cells_allowed_to_move<1
                number_cells_allowed_to_move=1;
            end
        end
    end
    disp(['Cells allowed to move: ',num2str(number_cells_allowed_to_move)]);
    disp(['Max redunding moves: ',num2str(max(max(history_map)))]);
    toc
end
disp('Converged !');
