% this is a file to look closer to the RB components that were created in
% createRB.m and compare truth/RB/DEIM solutions
clear all, clc

% units: geo or none
unit = 'geo';

%% paths/directories
% path to LaMEM executable
lamem = '"/home/chris/software/LaMEM/bin/opt/LaMEM"';
% add path to LaMEM matlab directory
addpath('/home/chris/software/LaMEM/matlab')
% path to src folder
addpath('/home/chris/Desktop/MA/RB_Stokes/reduced_basis_generation/src');
% LaMEM input file
input = '"Subduction2D.dat"';
% setup function file
setup = 'setup2D';
% name of folder to store RB components
folder = 'm6_30_40_max';


%% ======== length of the domain in x,y and z direction ===================
coordx = 700;
coordy = 20;
coordz = 700;


nel_x = 32;
nel_y = 2;
nel_z = 32;


n_p     =  nel_x*nel_y*nel_z;
n_vel   = ((nel_x+1)*nel_y*nel_z)+(nel_x*(nel_y+1)*nel_z)+(nel_x*nel_y*(nel_z+1));
n_tot = n_vel + n_p;

%% scaling factor to have velocities in cm/yr
length = 1000;
time   = 1e11;
cm_yr  = 0.01/(3600*24*365.25);

if (unit == 'geo')
    fac    = length/time/cm_yr;  % for geo units
else
    fac    = 1;                 % if nondimensional
end
%% load components
cd(folder);
load('B.mat');
load('eta_DEIM.mat');
load('M.mat');
load('M_DEIM.mat');
load('res_max.mat');
load('rho_DEIM.mat');
load('rhs_bl.mat');
load('rhs_bl_DEIM.mat');
cd ..

%% ================= check solutions ======================================

% example that lies in parameter space for RB problem
par1 = 32.456;
par2 = 36.342;

% markers
setup2D([par1,par2]);

% run simulation
[t1,t2] = system([lamem,' -ParamFile ', input]);
    
[t1,t2] = system([lamem,' -ParamFile ', input,' -only_matrix']);
%read data 
A         =  PetscBinaryRead('Matrices/Ass_A.bin');
Ml        =  PetscBinaryRead('Matrices/Ass_M.bin');
rhs       = -PetscBinaryRead('Matrices/rhs.bin');
rho       =  PetscBinaryRead('Matrices/rho.bin');
eta       =  PetscBinaryRead('Matrices/eta.bin');
sol_lamem =  PetscBinaryRead('Matrices/sol.bin');

% interpolate rho values
rho_i      = interpol_rho(nel_x, nel_y, nel_z, rho);
if (unit == 'geo')
rho_i = rho_i * 1e19;
end
rho_i_deim = rho_i;
rho_i      = nonzeros(rho_i);

%% truth solution from LaMEM
u_tl = sol_lamem;
u_tl(1:n_vel) = u_tl(1:n_vel)*fac;

%% truth solution with direct solver
J = A-Ml;
tic
u_t = J\rhs;
toc
u_t(1:n_vel) = u_t(1:n_vel)*fac;

%% RB solution without DEIM

% a few conts
m     = size(B(1,:),2);
eta   = [eta; 1];
N_K   = size(M(1,1,:),3);
K    = sparse(m,m);
N_R   = size(rhs_bl(1,:),2);
f    = sparse(m,1);

% assemble matrix with precomputed matrices
disp('direct solve with reduced basis by assembling the matrix with precomputed matrices:');
% assemble K
for i = 1:N_K
    K = K + (eta(i)*M(:,:,i));  
end

% assemble rhs
for i = 1:N_R
    f = f + (rho_i(i)*rhs_bl(:,i));
end

% solve RB system
tic
alpha = K\f;
u_RB = B * alpha;
toc
u_RB(1:n_vel) = u_RB(1:n_vel)*fac;

%% RB solution with DEIM
eta       =  PetscBinaryRead('Matrices/eta.bin');


% assemble matrix with precomputed  DEIM matrices
disp('direct solve with reduced basis with DEIM:');
m     = size(B(1,:),2);
N_K     = size(M_DEIM(1,1,:),3);
Kd      = sparse(m,m);
ct_eta  = eta_DEIM*eta;
ct_eta  = [ct_eta; 1];
N_R    = size(rhs_bl_DEIM(1,:),2);
fd     = sparse(m,1);
ct_rho = rho_DEIM*rho_i_deim;

for i = 1:N_K
    Kd = Kd + (ct_eta(i)*M_DEIM(:,:,i));  
end

for i = 1:N_R
    fd = fd + (ct_rho(i)*rhs_bl_DEIM(:,i));
end

tic
alphad = Kd\fd;
u_DEIM = B * alphad;
toc

u_DEIM(1:n_vel) = u_DEIM(1:n_vel)*fac;

%% calculate differences - overall
tr_RB  = abs(u_t)- abs(u_RB);
tr_D   = abs(u_t)- abs(u_DEIM);
RB_D   = abs(u_RB)- abs(u_DEIM);
max_trRB   = max(tr_RB);
max_tD    = max(tr_D);
max_RBD   = max(RB_D);

%% calculate differences - velocity
v_t  = u_t(1:n_vel);
v_RB = u_RB(1:n_vel);
v_D  = u_DEIM(1:n_vel);

v_tr_RB  = abs(v_t)- abs(v_RB);
v_tr_D   = abs(v_t)- abs(v_D);
v_RB_D   = abs(v_RB)- abs(v_D);
maxv_trRB   = max(v_tr_RB);
maxv_tD    = max(v_tr_D);
maxv_RBD   = max(v_RB_D);


%% ========= plot velocities =====================================================
% 
% plotting y velocity
figure(1)
% thruth solution
subplot(2,3,1)
sgtitle('z - velocity in xz plane');
[V3d_t] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, v_t,'z','xz');
x     = linspace(0,coordx,nel_x);
y     = linspace(0,coordz,nel_z+1);
[X,Y] = meshgrid(x,y);
pcolor(X,Y,V3d_t(:,:,2).'); c = colorbar;
c.Label.String = 'v_z [cm/yr]';
set(gca,'XTickLabel',[]);
shading interp;
title('thruth solution');
%xlabel('x');
ylabel('z [km]');
   
% RB velocity solution
subplot(2,3,2)
[V3d_RB] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, v_RB,'z','xz');
x     = linspace(0,coordx,nel_x);
y     = linspace(0,coordz,nel_z+1);
[X,Y] = meshgrid(x,y);
pcolor(X,Y,V3d_RB(:,:,2).'); c = colorbar;
c.Label.String = 'v_z [cm/yr]';
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
shading interp;
title('RB solution');
%xlabel('x');
%ylabel('z');

% DEIM solution
subplot(2,3,3)
[V3d_DEIM] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, v_D,'z','xz');
x     = linspace(0,coordx,nel_x);
y     = linspace(0,coordz,nel_z+1);
[X,Y] = meshgrid(x,y);
pcolor(X,Y,V3d_DEIM(:,:,2).'); c = colorbar;
c.Label.String = 'v_z [cm/yr]';
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
shading interp;
title('RB with DEIM');
%xlabel('x');
%ylabel('z');
    
% difference truth/RB
subplot(2,3,4)
diff_3D = V3d_t-V3d_RB;
%diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
[X,Y] = meshgrid(x,y);
pcolor(X,Y,diff_3D(:,:,2).'); c = colorbar;
c.Label.String = ' dv_z [cm/yr]';
shading interp;
title('difference btw thruth & RB');
xlabel('x [km]');
ylabel('z [km]');
    
% difference truth/DEIM
subplot(2,3,5)
diff_3D = V3d_t-V3d_DEIM;
%diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
[X,Y] = meshgrid(x,y);
pcolor(X,Y,diff_3D(:,:,2).'); c = colorbar;
c.Label.String = ' dv_z [cm/yr]';
set(gca,'YTickLabel',[]);
shading interp;
title('difference btw thruth & DEIM');
xlabel('x [km]');
%ylabel('z');
    
% difference RB/DEIM
subplot(2,3,6)
diff_3D = V3d_RB-V3d_DEIM;
%diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
[X,Y] = meshgrid(x,y);
pcolor(X,Y,diff_3D(:,:,2).'); c = colorbar;
c.Label.String = ' dv_z [cm/yr]';
set(gca,'YTickLabel',[]);
shading interp;
title('difference btw RB & DEIM');
xlabel('x [km]');
%ylabel('z');




%% plot residual   
%% plot max residual after adding a basis function
figure(2)
semilogy(res_max,'*-');
grid on;
ylabel('max residual');
xlabel('number of basis functions');

