clear all, clc
% add path to LaMEM matlab directory
addpath('/home/chris/software/LaMEM/matlab')
addpath('/home/chris/Desktop/MA/RB_Stokes')
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

%% ======== adjust RB parameters ==========================================
% viscosity of matrix
st   = 1;  % smallest parameter value
en   = 1;  % largest parameter value
n    = 1; % parameter spacing
par1 = linspace(st,en,n);

% viscosity of block
st   = 1;  % smallest parameter value
en   = 100; % largest parameter value
n    = 6;  % parameter spacing
par2 = linspace(st,en,n);


%% ======== reduced basis routine =========================================
[B,res, res_max, ETA] = create_RB(nel_x, nel_y, nel_z,par1,par2,1e-2);

%% === precompute matrixes that are needed to evaluate the RB solution ====
U       = [];
M       = precompute_mat(B, nel_x,nel_y,nel_z,0,ETA,U);

%% ======== DEIM ==========================================================
% apply DEIM to eta basis
ETA      = orth(ETA);      % orthonrmalize ETA matrix; important elsewise matrix tends to be singular
[U,P,p]  = DEIM(ETA);

M_DEIM   = precompute_mat(B, nel_x,nel_y,nel_z,1,U);
eta_DEIM = inv(P.'*U) * P.';


%% ================= check solutions ======================================
% create truth solution
eta1 = 1;
eta2 = 24;
system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(eta1),' -eta[1] ', num2str(eta2)]);
%[temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(loc1)),' -rho[1] ', num2str(par2(loc2))]);

% read data 
A         =  PetscBinaryRead('Matrices/Ass_A.bin');
M_lamem   =  PetscBinaryRead('Matrices/Ass_M.bin');
rhs       =  PetscBinaryRead('Matrices/rhs.bin');

% Solve linear system
[Sol,Sol_Vel,Sol_P,VV,VP,PV,PP,J] = solve_stokes(A,M_lamem,rhs,nel_x, nel_y, nel_z);

% solve RB
disp('direct solve of truth problem:');
tic 
u_truth = J\rhs;
toc

% assemble matrix with precomputed matrices
disp('direct solve with reduced basis by assembling the matrix with precomputed matrices:');
tic
m     = length(B(1,:));
eta   = PetscBinaryRead('Matrices/eta.bin');
N     = length(M(1,1,:));
K2    = sparse(m,m);
disp('Matrix assembling time:');
for i = 1:N
    K2 = K2 + (eta(i)*M(:,:,i));  
end
f2 = B.' * rhs;
alpha2 = K2\f2;
u_RB2 = B * alpha2;
toc

% assemble matrix with precomputed  DEIM matrices
disp('direct solve with reduced basis with DEIM:');

m     = length(B(1,:));
eta   =  PetscBinaryRead('Matrices/eta.bin');
N     = length(M_DEIM(1,1,:));
K3    = sparse(m,m);
tic
ct    = eta_DEIM*eta;
toc
tic
for i = 1:N
    K3 = K3 + (ct(i)*M_DEIM(:,:,i));  
end
toc
tic
f3 = B.' * rhs;
alpha3 = K3\f3;
u_DEIM = B * alpha3;
toc


%% calculate differences
ut      = u_truth(1:length(Sol_Vel));
urb     = u_RB(1:length(Sol_Vel));
urb2    = u_RB2(1:length(Sol_Vel));
uDEIM  = u_DEIM(1:length(Sol_Vel));
uRB_diff = max(ut-urb2)

%% ========= plot velocities =====================================================

% plotting y velocity
figure(1)
% thruth solution
subplot(2,3,1)
sgtitle('y - velocity in xz plane');
[V3d_t] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, ut,'y','xz');
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
[V3d_RB] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, urb2,'y','xz');
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
[V3d_DEIM] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, uDEIM,'y','xz');
x     = linspace(0,coordx,nel_x);
y     = linspace(0,coordz,nel_z);
[X,Y] = meshgrid(x,y);
pcolor(X,Y,V3d_DEIM(:,:,14).'); colorbar
shading interp;
title('DEIM solution');
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

    