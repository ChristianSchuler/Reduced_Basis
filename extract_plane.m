% input:
% nel_x    --> number of elements in x direction
% coord_x  --> length of the domain in x direction
% V_comp   --> x,y or z components of velocity
% plane    --> which plane to extract 
%  
% output:
% 3D structure of specific velocity component
function [V_3d] = extract_plane (nel_x, nel_y, nel_z, coordx, coordy, coordz, V_comp,plane)

%% y - z plane
if  plane == 'yz'

    V_3d = zeros(nel_y,nel_z,nel_x);
    
    for i = 1:nel_x 
        V_pl         = V_comp(i:nel_x:end);
        V_2d         = reshape(V_pl,nel_y,nel_z);   
        V_3d(:,:,i) = V_2d;   
    end
        
%% x - z plane        
elseif  plane == 'xz'   
    
    V_3d = zeros(nel_x,nel_z,nel_y);

    for i = 1:nel_y
        for j = 1:nel_x
            
             n_zplane      = nel_x * nel_y;
             V_col         = V_comp(1 + (i-1)*(nel_x)+(j-1):n_zplane:end);
             V_3d(j,:,i)   = V_col;
                                                        
            end
    end

%% x - z plane    
else  plane == 'xy'   
        

    V_3d = zeros(nel_x,nel_y,nel_z);

    for i = 1:nel_z
        n_zplane      = nel_x * nel_y;                                     
        V_pl          = V_comp(1+((i-1)*n_zplane):(n_zplane)+((i-1)*n_zplane));     
        V_2d          = reshape(V_pl,[nel_x,nel_y]);                           
        V_3d(:,:,i)   = V_2d;    
    end
%%     
end

end