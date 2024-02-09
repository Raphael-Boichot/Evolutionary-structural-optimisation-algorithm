function [boudary_conditions,growth,etching] = fun_ESO_algorithm(boudary_conditions,kp_k0,k0,heat_sink_temp,step_x,p_vol,rank_cells_to_exchange)
[height,width,~]=size(boudary_conditions);

%Listing the cells to grow
grow=zeros(height,width);
for k = 1:1:height
    for l = 1:1:width
        if boudary_conditions(k,l)==kp_k0
            if boudary_conditions(k+1,l)==k0; grow(k+1,l)=1; end
            if boudary_conditions(k-1,l)==k0; grow(k-1,l)=1; end
            if boudary_conditions(k,l+1)==k0; grow(k,l+1)=1; end
            if boudary_conditions(k,l-1)==k0; grow(k,l-1)=1; end
        end
    end
end

m=0;
for k = 1:1:height
    for l = 1:1:width
        if grow(k,l)==1
            m=m+1;
            grow_pos(m,1)=k;
            grow_pos(m,2)=l;
        end
    end
end

%Listing the cells to kill
etch=zeros(height,width);
for k = 1:1:height
    for l = 1:1:width
        if boudary_conditions(k,l)==k0
            if boudary_conditions(k+1,l)==kp_k0; etch(k+1,l)=1; end
            if boudary_conditions(k-1,l)==kp_k0; etch(k-1,l)=1; end
            if boudary_conditions(k,l+1)==kp_k0; etch(k,l+1)=1; end
            if boudary_conditions(k,l-1)==kp_k0; etch(k,l-1)=1; end
        end
    end
end

m=0;
for k = 1:1:height
    for l = 1:1:width
        if etch(k,l)==1
            m=m+1;
            etch_pos(m,1)=k;
            etch_pos(m,2)=l;
        end
    end
end

%identification of the cells to grow
boudary_conditions_temp=zeros(height,width,length(grow_pos));
for m=1:1:length(grow_pos)
    condi_temp=boudary_conditions;
    condi_temp(grow_pos(m,1),grow_pos(m,2))=kp_k0;
    boudary_conditions_temp(:,:,m)=condi_temp;
end
parfor m=1:length(grow_pos)
    % Variables output in this order :
    % 1. Distance of the hotest cell to the heat sink (scalar)
    % 2. Sum of cell entropy (scalar)
    % 3. Entropy map (matrix)
    % 4. Variance of temperatures accross the 1D adabatic borders (scalar)
    % 5. Variance of temperatures accross the 2D domain (scalar)
    % 6. Mean temperature (scalar)
    % 7. Maximal temperature accross the 2D domain (scalar)
    % 8. Map of temperatures (matrix)
    % 9. map of thermal gradients (matrix)
    % 10. Variance of gradients across the 2D domain (scalar)
    [~,~,~,~,~,~,temp_max,~,~,~]=finite_temp_direct_sparse(kp_k0,k0,heat_sink_temp,step_x,p_vol,boudary_conditions_temp(:,:,m));
    grow_pos(m,3)=temp_max;
end

grow_pos=sortrows(grow_pos,3);

%identification of the cells to kill
boudary_conditions_temp=zeros(height,width,length(etch_pos));
for m=1:1:length(etch_pos)
    condi_temp=boudary_conditions;
    condi_temp(etch_pos(m,1),etch_pos(m,2))=k0;
    boudary_conditions_temp(:,:,m)=condi_temp;
end
parfor m=1:length(etch_pos)
    [~,~,~,~,~,temp_max,~,~,~]=finite_temp_direct_sparse(kp_k0,k0,heat_sink_temp,step_x,p_vol,boudary_conditions_temp(:,:,m));
    etch_pos(m,3)=temp_max;
end
etch_pos=sortrows(etch_pos,3);

% growth=grow_pos(1:rank_cells_to_exchange,1:2);%list allowed to grow
% etching=etch_pos(1:rank_cells_to_exchange,1:2);%list allowed to be killed
% order_growth=randperm(number_cells_to_exchange);%shuffle data for more randomness
% order_etch=randperm(number_cells_to_exchange);%shuffle data for more randomness
% 
% %Cell exchange
% for i=1:1:number_cells_to_exchange
%     boudary_conditions(etching(order_etch(i),1),etching(order_etch(i),2))=k0;
%     boudary_conditions(growth(order_growth(i),1),growth(order_growth(i),2))=kp_k0;
% end

growth=grow_pos(ceil(rand*rank_cells_to_exchange),1:2);%list allowed to grow
etching=etch_pos(ceil(rand*rank_cells_to_exchange),1:2);%list allowed to be killed

boudary_conditions(etching(1,1),etching(1,2))=k0;
boudary_conditions(growth(1,1),growth(1,2))=kp_k0;