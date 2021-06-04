function [RHO_i] = interpol_rho (nel_x, nel_y, nel_z, RHO)

n_vx = ((nel_x+1)*nel_y*nel_z);
n_vy = (nel_x*(nel_y+1)*nel_z);
n_vz = (nel_x*nel_y*(nel_z+1));
n_p  = (nel_x*nel_y*nel_z);
n    = n_vx + n_vy + n_vz + n_p;
n_xy = nel_x*nel_y;
RHO_i = zeros(n,length(RHO(1,:)));
n_nz    = (nel_x*nel_y*(nel_z-1));

for j = 1:length(RHO(1,:))
    for i = 1:n_nz
        RHO_i((n_vx+n_vy+n_xy + i),j) = (RHO(i,j) + RHO(i+(n_xy),j));
    end
end

end