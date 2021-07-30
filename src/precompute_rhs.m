% function to compute rhs

function [rhs_final] = precompute_rhs (B, nel_x, nel_y, nel_z, g, mode, U, ip)

disp('=====================================================================')
disp('precompute rhs components ...');
disp('=====================================================================')

%% assemble rhs without rho values
% dimension of reduced basis
dB       = length(B(1,:));

n_vx = ((nel_x+1)*nel_y*nel_z);
n_vy = (nel_x*(nel_y+1)*nel_z);
n_vz = (nel_x*nel_y*(nel_z+1));
n_p  = (nel_x*nel_y*nel_z);
n    = n_vx + n_vy + n_vz + n_p;
n_xy = nel_x*nel_y;

m = n_vx-(2*n_xy);

% assemble matrix with affine decomposition of rhs (wo rho)
rhs_blank = sparse(n,m);
for i = 1:m
    rhs_blank((n_vx+n_vy+n_xy + i),i) = g/2;
end

B_t       =  B.';
B_t       = sparse(B_t);

if mode == 0
    
    rhs_final = sparse(dB,m);
    
    for i = 1:m
    rhs_final(:,i) = B_t * rhs_blank(:,i); 
    end

elseif mode == 1
    
    dim_DEIM = length(U(1,:));
    rhs_final = sparse(dB,dim_DEIM);
    
    for j = 1:dim_DEIM
        
        rhs_i = sparse(dB,1);
        
        for i = 1:m
        
            rhs_i = rhs_i + U(n_vx+n_vy+n_xy+i,j)* (B_t * rhs_blank(:,i)); 
        end
        
        rhs_final(:,j) = rhs_i;
        
    end
    
end
 
end