clear all, clc

% add path to LaMEM matlab directory
addpath('/home/chris/software/LaMEM/matlab')
addpath('/home/chris/Desktop/MA/RB_Stokes')

% number of elements
nel_x = 16;
nel_y = 16;
nel_z = 16;

% length of the domain in x,y and z direction
coordx = 1;
coordy = 1;
coordz = 1;

% create 'Matrices' folder or first delete it if it already exists
if not(isfolder('Matrices'))
    mkdir('Matrices');
else
    rmdir('Matrices','s');
    mkdir('Matrices');
end
addpath('/home/chris/Desktop/MA/RB_Stokes/reduced_basis_generation/Matrices')

%% ===================== adjust RB parameters =============================
% parameter space
% viscosity of matrix
st = 1;  % smallest parameter value
en = 1;  % largest parameter value
n  = 1; % parameter spacing
par1 = linspace(st,en,n);

% viscosity of block
st = 1;  % smallest parameter value
en = 1000; % largest parameter value
n  = 10;  % parameter spacing
par2 = linspace(st,en,n);
% =========================================================================

%% ==================== calculate RB decomposition matrices ===============
% collect all decompositionj matrices (usually many!!)
system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(1),' -eta[1] ', num2str(10),' -dump_decomposition_matrices']);


% assemble matrix
n_p     = nel_x*nel_y*nel_z;   % number of pressure nodes
A       =  PetscBinaryRead('Matrices/Ass_A.bin');
rhs     =  PetscBinaryRead('rhs.bin');
eta     =  PetscBinaryRead('eta.bin');
zrows   =  PetscBinaryRead('newzero.bin');
divV    = A(end-n_p+1:end,1:end-n_p);
divP    = A(1:end-n_p,end-n_p+1:end);

n_vis   = ((nel_x+1)*nel_y*nel_z)+(nel_x*(nel_y+1)*nel_z)+(nel_x*nel_y*(nel_z+1));
M       = sparse(n_vis,n_vis);
N       = n_p+((nel_x+1)*(nel_y+1)*nel_z)+((nel_x+1)*nel_y*(nel_z+1))+(nel_x*(nel_y+1)*(nel_z+1));

% tbd: 
for i = 1:N
    i_st    = num2str(i);
    Vis_fix = PetscBinaryRead(['Vis',i_st,'.bin']);
    Vis_fix = sparse(Vis_fix);
    Mi = eta(i)*Vis_fix;
    Mi = sparse(Mi);
    M  = M + Mi;
    i
end

% concatenate the different matrices
A_ass = zeros(n_vis+n_p,n_vis+n_p);% extract pressure and velocity divergence matrices from assembled A matrix
A_ass(1:n_vis,1:n_vis) = M;
A_ass(end-n_p+1:end,1:end-n_p) = divV;
A_ass(1:end-n_p,end-n_p+1:end) = divP;

% assemble boundary condition
for i = 1:length(zrows)
    i
    A_ass(zrows(i)+1,:)        = 0;
    A_ass(zrows(i)+1,zrows(i)+1) = 1;
end
% =========================================================================

%% reduced basis routine
[B,res, res_max] = create_RB(nel_x, nel_y, nel_z,par1,par2,1e-2);

%% DEIM for matrix
[U,P] = DEIM(B);

% precomputed DEIM part
 temp_inv = inv((P.'*U));   % inv((P.'*U))
 D = (U * temp_inv); %* P.');


%% create random truth solution and reduced basis solution

% create truth solution
eta1 = 1;
eta2 = 24;
[temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', eta1,' -eta[1] ', eta2]);
%[temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(loc1)),' -rho[1] ', num2str(par2(loc2))]);

% read data 
A   =  PetscBinaryRead('Mono_A.bin');
M   =  PetscBinaryRead('Mono_M.bin');
rhs =  PetscBinaryRead('r.1.bin');

% Solve linear system
[Sol,Sol_Vel,Sol_P,VV,VP,PV,PP,J] = solve_stokes(A,M,rhs,nel_x, nel_y, nel_z);

% solve RB
disp('direct solve of truth problem:');
tic 
Sol_test = J\rhs;
toc

disp('direct solve with reduced basis:');
tic
K = B.' * (J) * B;
f = B.' * rhs;
alpha = K\f;
u_RB = B * alpha;
toc

% disp('greedy - DEIM:');
% tic
% K = (B.' * D) * (P.' * (A * B));
% f = B.' * rhs;
% %K = B.' * (A) * B;
% %f = B.' * rhs;
% alpha = K\f;
% u_RB = B * alpha;
% toc


%% plotting

% plotting y velocity
    figure(1)
    % truth velocity solution
    subplot(2,3,1)
    [V_3d_t] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, Sol_Vel,'y','xz');
    x     = linspace(0,coordx,nel_x);
    y     = linspace(0,coordz,nel_z);
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,V_3d_t(:,:,14).'); colorbar
    shading interp;
    title('y - velocity in xz plane (truth solution)');
    xlabel('x');
    ylabel('z');
    
    % RB velocity solution
    subplot(2,3,2)
    [V_3d_RB] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, u_RB(1:length(Sol_Vel)),'y','xz');
    x     = linspace(0,coordx,nel_x);
    y     = linspace(0,coordz,nel_z);
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,V_3d_RB(:,:,14).'); colorbar
    shading interp;
    title('y - velocity in xz plane (RB solution)');
    xlabel('x');
    ylabel('z');
    
     % difference
    subplot(2,3,3)
    diff_3D = V_3d_t-V_3d_RB;
    diff_3D = diff_3D;%/max(max(V_3d_t(:,:,14)));
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,diff_3D(:,:,14).'); colorbar
    shading interp;
    title('difference between both solutions (normalized)');
    xlabel('x');
    ylabel('z');
    
      subplot(2,3,4)
    [V_3d_t] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, Sol_Vel,'z','xz');
    x     = linspace(0,coordx,nel_x);
    y     = linspace(0,coordz,nel_z+1);
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,V_3d_t(:,:,14).'); colorbar
    shading interp;
    title('z - velocity in xz plane (truth solution)');
    xlabel('x');
    ylabel('z');
    
    % RB velocity solution
    subplot(2,3,5)
    [V_3d_RB] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, u_RB(1:length(Sol_Vel)),'z','xz');
    x     = linspace(0,coordx,nel_x);
    y     = linspace(0,coordz,nel_z+1);
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,V_3d_RB(:,:,14).'); colorbar
    shading interp;
    title('z - velocity in xz plane (RB solution)');
    xlabel('x');
    ylabel('z');
    
     % difference
    subplot(2,3,6)
    diff_3D = V_3d_t-V_3d_RB;
    diff_3D = diff_3D;%/max(max(V_3d_t(:,:,14)));
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,diff_3D(:,:,14).'); colorbar
    shading interp;
    title('difference between both solutions (normalized)');
    xlabel('x');
    ylabel('z');
    
    
%     % truth pressure solution
%     subplot(2,3,4)
%     [P_3d_t] = extract_plane (nel_x, nel_y, nel_z, coordx, coordy, coordz, Sol_P,'xy');
%     x     = linspace(0,coordx,nel_x);
%     y     = linspace(0,coordy,nel_y);
%     [X,Y] = meshgrid(x,y);
%     pcolor(X,Y,P_3d_t(:,:,14).'); colorbar
%     shading interp;
%     title('pressure ijn x-z plane');
%     xlabel('x');
%     ylabel('z');
%     
%     % RB velocity solution
%     subplot(2,3,5)
%     [P_3d_RB] = extract_plane (nel_x, nel_y, nel_z, coordx, coordy, coordz, u_RB(length(Sol_Vel)+1:end),'xy');
%     x     = linspace(0,coordx,nel_x);
%     y     = linspace(0,coordz,nel_z);
%     [X,Y] = meshgrid(x,y);
%     pcolor(X,Y,P_3d_RB(:,:,14).'); colorbar
%     shading interp;
%     title('pressure ijn x-z plane');
%     xlabel('x');
%     ylabel('z');
%     
%      % difference
%     subplot(2,3,6)
%     diff_3D = P_3d_t-P_3d_RB;
%     [X,Y] = meshgrid(x,y);
%     pcolor(X,Y,diff_3D(:,:,14).'); colorbar
%     shading interp;
%     title('difference between both solutions');
%     xlabel('x');
%     ylabel('z');
    
    
% plot max residual after adding a basis function
figure(2)
semilogy(res_max,'*-');
grid on;
ylabel('max residual');
xlabel('number of basis functions');


% max errors
max_vx = max(abs(Sol_Vel(1:4352) - u_RB(1:4352)));
max_vy = max(abs(Sol_Vel(4352:4352*2) - u_RB(4352:4352*2)));
max_vz = max(abs(Sol_Vel(4352*2:end) - u_RB(4352*2:4352*3)));
max_P = max(abs(Sol_P - u_RB(length(Sol_Vel)+1:end)));
    
