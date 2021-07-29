clear all, clc

%% paths/directories
% path to external packages
addpath(genpath('/home/chris/Desktop/MA/RB_Stokes/reduced_basis_generation/external_packages'));
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

% create 'Matrices' folder or first delete it if it already exists
if not(isfolder('Matrices'))
    mkdir('Matrices');
else
    rmdir('Matrices','s');
    mkdir('Matrices');
end

%% ======== number of elements ============================================
nel_x  = 32;
nel_y  = 2;
nel_z  = 32;

%% ======== length of the domain in x,y and z direction ===================
coordx = 700;
coordy = 20;
coordz = 700;

%% ======= graviational accleration =======================================
g    = -9.81;   % acceleration in z direction; assuming that acceleration of other directions is 0

%% ======== adjust RB parameters ==========================================
% define parameters/their ranges and distribution
% parameter 1
st   = 20;  % smallest parameter value
en   = 60;  % largest parameter value
n    = 8; % parameter spacing
par1 = linspace(st,en,n);
% parameter 2
st   = 20;  % smallest parameter value
en   = 60;  % largest parameter value
n    = 8; % parameter spacing
par2 = linspace(st,en,n);

% allcomb creates all possible combinations of input parameters
par = allcomb(par1,par2);

%% ======== reduced basis routine =========================================
tol  = 1e-2;    % tolerance for Greedy algorithm
n    = 10;      % orthogonalization of basis every n steps
[B, res_max, ETA, RHO] = Reduced_Basis(lamem, input,setup, nel_x, nel_y, nel_z,par,tol,n);

% %% ======== extract decomposition matrices ================================
% preMat = extract_preMat(lamem, input, nel_x,nel_y,nel_z);
% preMat(:,:,end) = [];  % last matrix was only built for scaling --> dump it
 
%% decomposition matrices without DEIM
% U       = [];
% ip      = [];
% M      = precompute_mat(lamem, input2, preMat, B, nel_x, nel_y, nel_z, 0, U);
% rhs_bl  = precompute_rhs(B, nel_x, nel_y, nel_z, g , 0, U,ip);

%{
%% ======== precomputed Jacobian matrixes with DEIM =======================
% apply DEIM to eta basis
ETA      = orth(ETA);      % orthonormalize ETA matrix; important elsewise matrix tends to be singular
[U,P,p]  = DEIM(ETA);
eta_DEIM = inv(P.'*U) * P.';
M_DEIM   = precompute_mat(lamem, input2, preMat, B, nel_x, nel_y, nel_z, 1, U);

RHO_i    = interpol_rho(nel_x, nel_y, nel_z, RHO);

RHO_i      = orth(RHO_i);  % orthonormalize ETA matrix; important elsewise matrix tends to be singular
[U,P,p]    = DEIM(RHO_i);
[ip, col] = find(P);
rho_DEIM   = inv(P.'*U) * P.';

rhs_bl_DEIM = precompute_rhs(B, nel_x, nel_y, nel_z, g , 1, U, ip);

disp('reduced basis offline computations took:'); 
toc
% 
% % create 'Matrices' folder or first delete it if it already exists
% if not(isfolder('RB_components'))
%     mkdir('RB_components');
% else
%     rmdir('RB_components','s');
%     mkdir('RB_components');
% end
% 
% save('saveM_DEIM.mat','M_DEIM');
%}
%% ================= check solutions ======================================
%create truth solution
% eta = 66;
% rho = 10;
% [t1, t2] = system([lamem,' -ParamFile ../', input, ' -eta[1] ', num2str(eta),' -rho[1] ', num2str(rho)]);
% 
% inputG      = input;
% inputG(1)   = [];
% inputG(end) = [];


% copyfile(inputG,'geometry.dat');
% radius 		= 	0.12;	
% str = ['<SphereStart>' newline 'phase  = 1' newline 'center = 0.5 0.5 0.5' newline 'radius = ' num2str(radius) newline '<SphereEnd>'];
% fid=fopen('geometry.dat','a+');
% fprintf(fid, str);
% fclose(fid);
% [t1, t2] = system([lamem,' -ParamFile geometry.dat']);

% setup3D(194,127);
% % run simulation
% system([lamem,' -ParamFile Subduction3D.dat']);
%     
% system([lamem,' -ParamFile Subduction3D.dat -only_matrix']);
% %read data 
% A         =  PetscBinaryRead('Matrices/Ass_A.bin');
% M_lamem   =  PetscBinaryRead('Matrices/Ass_M.bin');
% rhs       =  -PetscBinaryRead('Matrices/rhs.bin');
% rho       =  PetscBinaryRead('Matrices/rho.bin');
% sol_lamem =  PetscBinaryRead('Matrices/sol.bin');
% 
% J = A-M_lamem;
% u_sol = J\rhs;
%    
% n_nz    = (nel_x*nel_y*(nel_z-1));
% n_vy = (nel_x*(nel_y+1)*nel_z);
% n_vx = ((nel_x+1)*nel_y*nel_z);
% n_xy    = nel_x*nel_y;
% n_p     =  nel_x*nel_y*nel_z;
% rho_i = zeros(n_nz,1);
% 
%  % total number of velocity nodes in the mesh
% n_velx = (nel_x+1)*nel_y*nel_z;
% n_vely = (nel_y+1)*nel_x*nel_z;
% n_velz = (nel_z+1)*nel_x*nel_y;
% n_vel  = n_velx + n_vely + n_velz; % total number of all velocity points
% 
% % pressure nodes
% n_p = nel_x * nel_y * nel_z;
% 
% %total nodes
% N = n_vel+n_p;
% 
% for i = 1:n_nz
%     rho_i(i) = (rho(i) + rho(i+(n_xy)));
% end
% 
% rho_i_deim = zeros(N,1);
% 
% for i = 1:n_nz
%     rho_i_deim((n_vx+n_vy+n_xy + i)) = (rho(i) + rho(i+(n_xy)));
% end
% 
% % rhs wo rho
% g_fac = g/2;
% 
% % calculate rhs
%             rho_i = zeros(n_nz,1);
% 
%             for i = 1:n_nz
%                 rho_i(i) = (rho(i) + rho(i+(n_xy)));
%             end
%             
%             rhs = sparse(N,1);
% 
%             for i = 1:(n_velz-(2*n_xy))
%                 rhs(n_velx+n_vely+n_xy+i) = rho_i(i)*g_fac;
%             end
% 
% % % Solve linear system
% %[Sol,Sol_Vel,Sol_P,VV,VP,PV,PP,J] = solve_stokes(A,M_lamem,rhs,nel_x, nel_y, nel_z);
% 
% % solve RB
% disp('direct solve of truth problem:');
% tic 
% u_truth1 = J\(rhs);
% toc
% 
% disp('direct solve with reduced basis:');
% tic
% K = B.' * (J) * B;
% f = B.' * rhs;
% alpha = K\f;
% u_RB = B * alpha;
% toc
% 
% eta   = PetscBinaryRead('Matrices/eta.bin');
% % assemble matrix with precomputed matrices
% disp('direct solve with reduced basis by assembling the matrix with precomputed matrices:');
% tic
% m     = length(B(1,:));
% eta   = [eta; 1; 1];
% N_K     = length(M(1,1,:));
% K2    = sparse(m,m);
% disp('Matrix assembling time:');
% % assemble K
% for i = 1:N_K
%     K2 = K2 + (eta(i)*M(:,:,i));  
% end
% 
% % assemble rhs
% N_R   = length(rhs_bl(1,:));
% f2    = sparse(m,1);
% for i = 1:N_R
%     f2 = f2 + (rho_i(i)*rhs_bl(:,i));
% end
% 
% alpha2 = K2\f2;
% u_RB2 = B * alpha2;
% toc
% 
% eta   =  PetscBinaryRead('Matrices/eta.bin');
% % assemble matrix with precomputed  DEIM matrices
% disp('direct solve with reduced basis with DEIM:');
% tic
% m     = length(B(1,:));
% % assemble matrix K
% N_K     = length(M_DEIM(1,1,:));
% K3      = sparse(m,m);
% ct_eta  = eta_DEIM*eta;
% ct_eta  = [ct_eta; 1; 1];
% for i = 1:N_K
%     K3 = K3 + (ct_eta(i)*M_DEIM(:,:,i));  
% end
% 
% % assemble rhs
% N_R    = length(rhs_bl_DEIM(1,:));
% f3     = sparse(m,1);
% ct_rho = rho_DEIM*rho_i_deim;
% for i = 1:N_R
%     f3 = f3 + (ct_rho(i)*rhs_bl_DEIM(:,i));
% end
% 
% alpha3 = K3\f3;
% u_DEIM = B * alpha3;
% toc
% 

% %% calculate differences
% %ut       = u_truth1(1:length(Sol_Vel));
% u_lamem  = sol_lamem(1:n_vel);
% urb      = u_RB(1:n_vel);
% urb2     = u_RB2(1:length(Sol_Vel));
% uDEIM    = u_DEIM(1:length(Sol_Vel));
% % uRB_diff = max(ut-urb2);

%rmdir('Matrices','s');

%% ========= plot velocities =====================================================
% 
% % plotting y velocity
% figure(1)
% % thruth solution
% subplot(2,3,1)
% sgtitle('z - velocity in xz plane');
% [V3d_t] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, u_lamem,'z','yz');
% x     = linspace(0,coordx,nel_x)*800;
% y     = linspace(0,coordz,nel_z+1)*690;
% [X,Y] = meshgrid(x,y);
% pcolor(X,Y,V3d_t(:,:,14).'); colorbar
% 
% 
% shading interp;
% title('thruth solution');
% xlabel('x');
% ylabel('z');
%    
% % RB velocity solution
% subplot(2,3,2)
% [V3d_RB] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, urb,'z','yz');
% x     = linspace(0,coordx,nel_x)*800;
% y     = linspace(0,coordz,nel_z+1)*690;
% [X,Y] = meshgrid(x,y);
% pcolor(X,Y,V3d_RB(:,:,14).'); colorbar
% 
% 
% shading interp;
% title('RB solution');
% xlabel('x');
% ylabel('z');
% 
% % difference truth/RB
% subplot(2,3,3)
% diff_3D = V3d_t-V3d_RB;
% %diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
% [X,Y] = meshgrid(x,y);
% pcolor(X,Y,diff_3D(:,:,14).'); colorbar
% 
% shading interp;
% title('difference btw thruth & RB');
% xlabel('x');
% ylabel('z');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % DEIM solution
% subplot(2,3,3)
% [V3d_DEIM] = arrange_vel (nel_x, nel_y, nel_z, coordx, coordy, coordz, uDEIM,'z','xy');
% x     = linspace(0,coordx,nel_x);
% y     = linspace(0,coordz,nel_z);
% [X,Y] = meshgrid(x,y);
% pcolor(X,Y,V3d_DEIM(:,:,14).'); colorbar
% shading interp;
% title('RB with DEIM');
% xlabel('x');
% ylabel('z');
%     
% % difference truth/RB
% subplot(2,3,4)
% diff_3D = V3d_t-V3d_RB;
% %diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
% [X,Y] = meshgrid(x,y);
% pcolor(X,Y,diff_3D(:,:,14).'); colorbar
% shading interp;
% title('difference btw thruth & RB');
% xlabel('x');
% ylabel('z');
%     
% % difference truth/DEIM
% subplot(2,3,5)
% diff_3D = V3d_t-V3d_DEIM;
% %diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
% [X,Y] = meshgrid(x,y);
% pcolor(X,Y,diff_3D(:,:,14).'); colorbar
% shading interp;
% title('difference btw thruth & DEIM');
% xlabel('x');
% ylabel('z');
%     
% % difference RB/DEIM
% subplot(2,3,6)
% diff_3D = V3d_RB-V3d_DEIM;
% %diff_3D = diff_3D/max(max(V_3d_t(:,:,14)));
% [X,Y] = meshgrid(x,y);
% pcolor(X,Y,diff_3D(:,:,14).'); colorbar
% shading interp;
% title('difference btw RB & DEIM');
% xlabel('x');
% ylabel('z');
% 
% %}


%% plot residual   
%% plot max residual after adding a basis function
figure(2)
semilogy(res_max,'*-');
grid on;
ylabel('max residual');
xlabel('number of basis functions');












