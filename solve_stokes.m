function [Sol,Sol_Vel,Sol_P,VV,VP,PV,PP,J] = solve_stokes(A,M,rhs,nel_x, nel_y, nel_z)


% Jacobian
J = A - M;

% total number of velocity nodes in the mesh
n_velx = (nel_x+1)*nel_y*nel_z;
n_vely = (nel_y+1)*nel_x*nel_z;
n_velz = (nel_z+1)*nel_x*nel_y;
n_vel  = n_velx + n_vely + n_velz; % total number of all velocity points

% collect block matrixes
VV = J(1:n_vel,1:n_vel);
VP = J(n_vel+1:end,1:n_vel);
PV = J(1:n_vel,n_vel+1:end);
PP = J(n_vel+1:end,n_vel+1:end);

% pressure nodes
n_p = nel_x * nel_y * nel_z;

% Solve with a direct solver
Sol = J\rhs;                 % Solve
Sol_Vel         =   Sol(1:n_vel);
Sol_P           =   Sol(n_vel + 1:end);

end