%https://github.com/Raphael-Boichot/Evolutionary-structural-optimisation-algorithm
function initial_boundary_conditions=init_image(initial_boundary_conditions,conductive_cells, k0, kp_k0)
[height,width,~]=size(initial_boundary_conditions);
k=0;

while k<conductive_cells
    h=ceil(rand*height);
    l=ceil(rand*width);
    
    if initial_boundary_conditions(h,l)==k0 
    initial_boundary_conditions(h,l)=k0*kp_k0; 
    k=k+1;
    end
        
end