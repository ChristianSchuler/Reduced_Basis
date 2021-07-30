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
g    = -9.81;   % acceleration in z direction;acceleration of other directions -> 0

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
tol  = 1e-3;    % tolerance for Greedy algorithm
n    = 10;      % orthogonalization of basis every n steps
[B, res_max, ETA, RHO] = Reduced_Basis(lamem, input,setup, nel_x, nel_y, nel_z,par,tol,n);

%% ======== extract decomposition matrices ================================
preMat = extract_preMat(lamem, input, nel_x,nel_y,nel_z);
preMat(:,:,end) = [];  % last matrix was only built for scaling --> dump it
 
%% decomposition matrices without DEIM
U       = [];
ip      = [];
M       = precompute_mat(preMat, B, nel_x, nel_y, nel_z, 0, U);
rhs_bl  = precompute_rhs(B, nel_x, nel_y, nel_z, g , 0, U,ip);

% %{
%% ======== precomputed Jacobian matrixes with DEIM =======================
% apply DEIM to eta basis
ETA      = orth(ETA);      % orthonormalize ETA matrix; important elsewise matrix tends to be singular
[U,P,p]  = DEIM(ETA);
eta_DEIM = inv(P.'*U) * P.';
M_DEIM   = precompute_mat(preMat, B, nel_x, nel_y, nel_z, 1, U);

RHO_i    = interpol_rho(nel_x, nel_y, nel_z, RHO);

RHO_i      = orth(RHO_i);  % orthonormalize ETA matrix; important elsewise matrix tends to be singular
[U,P,p]    = DEIM(RHO_i);
[ip, col] = find(P);
rho_DEIM   = inv(P.'*U) * P.';

rhs_bl_DEIM = precompute_rhs(B, nel_x, nel_y, nel_z, g , 1, U, ip);

disp('reduced basis offline computations took:'); 
toc

% create 'Matrices' folder or first delete it if it already exists
if not(isfolder('RB_components'))
    mkdir('RB_components');
else
    rmdir('RB_components','s');
    mkdir('RB_components');
end

save('saveM_DEIM.mat','M_DEIM');
%}











