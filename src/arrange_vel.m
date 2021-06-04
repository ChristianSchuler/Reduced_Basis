% input:
% nel_x    --> number of elements in x direction
% coord_x  --> length of the domain in x direction
% Sol_V    --> velocity solution vector
% vel_comp --> which velocity component shoul be plotted
% plane    --> which plane should be plotted
%  
% output:
% 3d structure of specific velocity component
function [V_3d] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, Sol_Vel,vel_comp,plane)


%% x velocity

if     vel_comp == 'x'
     
    % extract x velocities
    V_x = Sol_Vel(1:((nel_x+1)*nel_y*nel_z)); 
    
    % one extra velocity component in x direction
    nel_x = nel_x + 1; 
    
    if      plane == 'yz'
        V_3d = extract_plane(nel_x, nel_y, nel_z, coordx, coordy, coordz, V_x,'yz');  
        
    elseif  plane == 'xz'  
        V_3d = extract_plane(nel_x, nel_y, nel_z, coordx, coordy, coordz, V_x,'xz');  
    
    else    plane == 'xy'   
        V_3d = extract_plane(nel_x, nel_y, nel_z, coordx, coordy, coordz, V_x,'xy');  
        
    end 
 
%% y velocity    
elseif vel_comp == 'y'
 
    % extract y velocities
    sV_x = (nel_x+1)*nel_y*nel_z; % size V_x
    V_y = Sol_Vel((sV_x+1):sV_x+(nel_x*(nel_y+1)*nel_z));
    
    % one extra velocity component in y direction
    nel_y = nel_y + 1;
    
    if      plane == 'yz'
        V_3d = extract_plane(nel_x, nel_y, nel_z, coordx, coordy, coordz, V_y,'yz');
        
    elseif  plane == 'xz'   
        V_3d = extract_plane(nel_x, nel_y, nel_z, coordx, coordy, coordz, V_y,'xz'); 
    
    else    plane == 'xy'   
        V_3d = extract_plane(nel_x, nel_y, nel_z, coordx, coordy, coordz, V_y,'xy'); 
    end
        
%% z velocity       
else  vel_comp == 'z'
   
    % extract z velocities
    sV_x = (nel_x+1)*nel_y*nel_z; % size V_x
    sV_y = nel_x*(nel_y+1)*nel_z; % size V_y
    V_z = Sol_Vel((sV_x+sV_y+1):end);
    
    % one extra velocity component in y direction
    nel_z = nel_z + 1;
    
    if      plane == 'yz'
        V_3d = extract_plane(nel_x, nel_y, nel_z, coordx, coordy, coordz, V_z,'yz');
        
    elseif  plane == 'xz'   
        V_3d = extract_plane(nel_x, nel_y, nel_z, coordx, coordy, coordz, V_z,'xz'); 
    
    else    plane == 'xy'   
        V_3d = extract_plane(nel_x, nel_y, nel_z, coordx, coordy, coordz, V_z,'xy'); 
        
    end
%%        
end
    
end