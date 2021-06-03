clear all, clc
% add path to LaMEM matlab directory
addpath('/home/chris/software/LaMEM/matlab')
addpath('/home/chris/Desktop/MA/RB_Stokes')
addpath('/home/chris/Desktop/MA/RB_Stokes/reduced_basis_generation/ndSparse');

% create 'Matrices' folder or first delete it if it already exists
if not(isfolder('Matrices'))
    mkdir('Matrices');
else
    rmdir('Matrices','s');
    mkdir('Matrices');
end

%% ======== number of elements ============================================
nel_x  = 16;
nel_y  = 16;
nel_z  = 16;

%% ======== length of the domain in x,y and z direction ===================
coordx = 1;
coordy = 1;
coordz = 1;

%% ======= graviational accleration =======================================
g = -1; % acceleration in z direction; assuming that acceleration of other directions is 0

%% ======== adjust RB parameters ==========================================
% viscosity of block
st   = 10;  % smallest parameter value
en   = 100;  % largest parameter value
n    = 4; % parameter spacing
par1 = linspace(st,en,n);

% density of block
st   = 10;  % smallest parameter value
en   = 10; % largest parameter value
n    = 1;  % parameter spacing
par2 = linspace(st,en,n);


%% ======== reduced basis routine =========================================
[B, res_max, ETA, RHO] = create_RB(nel_x, nel_y, nel_z,g,par1,par2,1e-3);
%B = orth(B);

%% extract decomposition matrices
preMat = extract_preMat (nel_x,nel_y,nel_z);

%% ======== precompute matrixes from Jacobian =============================
U       = [];
M       = precompute_mat(preMat, B, nel_x, nel_y, nel_z, 0, U);

%% ======== precomputed Jacobian matrixes with DEIM =======================
% apply DEIM to eta basis
ETA      = orth(ETA);      % orthonormalize ETA matrix; important elsewise matrix tends to be singular
[U,P,p]  = DEIM(ETA);
eta_DEIM = inv(P.'*U) * P.';

M_DEIM   = precompute_mat(preMat, B, nel_x, nel_y, nel_z, 1, U);


%% ======= DEIM with rhs ==================================================
% this method holds for rhs where only the gravitational potential matters 
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

RHO_i    = orth(RHO_i);      % orthonormalize ETA matrix; important elsewise matrix tends to be singular
[U,P,p]  = DEIM(RHO_i);
rho_DEIM = inv(P.'*U) * P.';

[rhs_bl]      = precompute_rhs(B, nel_x, nel_y, nel_z, g , 0, U);
[rhs_bl_DEIM] = precompute_rhs(B, nel_x, nel_y, nel_z, g , 1, U);

%% ================= check solutions ======================================
% create truth solution
eta = 66;
rho = 10;
system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect2.dat -eta[1] ', num2str(eta),' -rho[1] ', num2str(rho)]);

% read data 
A         =  PetscBinaryRead('Matrices/Ass_A.bin');
M_lamem   =  PetscBinaryRead('Matrices/Ass_M.bin');
rho       =  PetscBinaryRead('Matrices/rho.bin');
sol_lamem =  PetscBinaryRead('Matrices/sol.bin');
   
n_nz    = (nel_x*nel_y*(nel_z-1));
n_xy    = nel_x*nel_y;
rho_i = zeros(n_nz,1);

for i = 1:n_nz
    rho_i(i) = (rho(i) + rho(i+(n_xy)));
end

rho_i_deim = zeros(n,1);

for i = 1:n_nz
    rho_i_deim((n_vx+n_vy+n_xy + i)) = (rho(i) + rho(i+(n_xy)));
end

% rhs wo rho
g_fac = g/2;

% variables for calculating rhs
n_nz    = (nel_x*nel_y*(nel_z-1));
n_xy    = nel_x*nel_y;

 % total number of velocity nodes in the mesh
n_velx = (nel_x+1)*nel_y*nel_z;
n_vely = (nel_y+1)*nel_x*nel_z;
n_velz = (nel_z+1)*nel_x*nel_y;
n_vel  = n_velx + n_vely + n_velz; % total number of all velocity points

% pressure nodes
n_p = nel_x * nel_y * nel_z;

%total nodes
N = n_vel+n_p;
% calculate rhs
            rho_i = zeros(n_nz,1);

            for i = 1:n_nz
                rho_i(i) = (rho(i) + rho(i+(n_xy)));
            end
            
            rhs = sparse(N,1);

            for i = 1:n_velz-(2*n_xy)
                rhs(n_velx+n_vely+n_xy+i) = rho_i(i)*g_fac;
            end

% % Solve linear system
[Sol,Sol_Vel,Sol_P,VV,VP,PV,PP,J] = solve_stokes(A,M_lamem,rhs,nel_x, nel_y, nel_z);

% solve RB
disp('direct solve of truth problem:');
tic 
u_truth1 = J\(-rhs);
toc

disp('direct solve with reduced basis:');
tic
K = B.' * (J) * B;
f = B.' * rhs;
alpha = K\f;
u_RB = B * alpha;
toc

% assemble matrix with precomputed matrices
disp('direct solve with reduced basis by assembling the matrix with precomputed matrices:');
tic
m     = length(B(1,:));
eta   = PetscBinaryRead('Matrices/eta.bin');
eta   = [eta; 1; 1];
N_K     = length(M(1,1,:));
K2    = sparse(m,m);
disp('Matrix assembling time:');
% assemble K
for i = 1:N_K
    K2 = K2 + (eta(i)*M(:,:,i));  
end

% assemble rhs
N_R   = length(rhs_bl(1,:));
f2    = sparse(m,1);
for i = 1:N_R
    f2 = f2 + (rho_i(i)*rhs_bl(:,i));
end

%f2 = B.' * -rhs;
alpha2 = K2\(f2);
u_RB2 = B * alpha2;
toc

% assemble matrix with precomputed  DEIM matrices
disp('direct solve with reduced basis with DEIM:');
tic
m     = length(B(1,:));
eta   =  PetscBinaryRead('Matrices/eta.bin');
% assemble matrix K
N_K     = length(M_DEIM(1,1,:));
K3      = sparse(m,m);
ct_eta  = eta_DEIM*eta;
ct_eta  = [ct_eta; 1; 1];
for i = 1:N_K
    K3 = K3 + (ct_eta(i)*M_DEIM(:,:,i));  
end



% % assemble rhs
% N_R    = length(rhs_bl_DEIM(1,:));
% f3     = sparse(m,1);
% ct_rho = rho_DEIM*rho_i_deim;
% for i = 1:N_R
%     f3 = f3 + (ct_rho(i)*rhs_bl_DEIM(:,i));
% end

f3 = B.' * rhs;
alpha3 = K3\f3;
u_DEIM = B * alpha3;
toc


%% calculate differences
%ut       = u_truth1(1:length(Sol_Vel));
u_lamem  = sol_lamem(1:length(Sol_Vel));
urb      = u_RB(1:length(Sol_Vel));
urb2     = u_RB2(1:length(Sol_Vel));
uDEIM    = u_DEIM(1:length(Sol_Vel));
% uRB_diff = max(ut-urb2);

%% ========= plot velocities =====================================================

% plotting y velocity
figure(1)
% thruth solution
subplot(2,3,1)
sgtitle('y - velocity in xz plane');
[V3d_t] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, u_lamem,'z','xy');
x     = linspace(0,coordx,nel_x);
y     = linspace(0,coordz,nel_z);
[X,Y] = meshgrid(x,y);
pcolor(X,Y,V3d_t(:,:,14).'); colorbar

shading interp;
title('thruth solution');
xlabel('x');
ylabel('z');
   
% RB velocity solution
subplot(2,3,2)
[V3d_RB] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, urb2,'z','xy');
x     = linspace(0,coordx,nel_x);
y     = linspace(0,coordz,nel_z);
[X,Y] = meshgrid(x,y);
pcolor(X,Y,V3d_RB(:,:,14).'); colorbar
shading interp;
title('RB solution');
xlabel('x');
ylabel('z');
    
% DEIM solution
subplot(2,3,3)
[V3d_DEIM] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, uDEIM,'z','xy');
x     = linspace(0,coordx,nel_x);
y     = linspace(0,coordz,nel_z);
[X,Y] = meshgrid(x,y);
pcolor(X,Y,V3d_DEIM(:,:,14).'); colorbar
shading interp;
title('RB with DEIM');
xlabel('x');
ylabel('z');
    
% difference truth/RB
subplot(2,3,4)
diff_3D = V3d_t-V3d_RB;
%diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
[X,Y] = meshgrid(x,y);
pcolor(X,Y,diff_3D(:,:,14).'); colorbar
shading interp;
title('difference btw thruth & RB');
xlabel('x');
ylabel('z');
    
% difference truth/DEIM
subplot(2,3,5)
diff_3D = V3d_t-V3d_DEIM;
%diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
[X,Y] = meshgrid(x,y);
pcolor(X,Y,diff_3D(:,:,14).'); colorbar
shading interp;
title('difference btw thruth & DEIM');
xlabel('x');
ylabel('z');
    
% difference RB/DEIM
subplot(2,3,6)
diff_3D = V3d_RB-V3d_DEIM;
%diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
[X,Y] = meshgrid(x,y);
pcolor(X,Y,diff_3D(:,:,14).'); colorbar
shading interp;
title('difference btw RB & DEIM');
xlabel('x');
ylabel('z');

%% plot residual   
% plot max residual after adding a basis function
figure(2)
semilogy(res_max,'*-');
grid on;
ylabel('max residual');
xlabel('number of basis functions');

    