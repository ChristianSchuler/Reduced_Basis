% script to compare assembled and dumped matrix
clear all,clc

format long;

%% paths/directories
% add path to LaMEM matlab directory
addpath('/home/chris/software/LaMEM/matlab')
% path to src folder
addpath('/home/chris/Desktop/MA/RB_Stokes/reduced_basis_generation/src');
% path to ndSparse package
addpath('/home/chris/Desktop/MA/RB_Stokes/reduced_basis_generation/ndSparse');
% path to mtimessx package
addpath('/home/chris/Desktop/MA/RB_Stokes/reduced_basis_generation/mtimesx_20110223');
% path to LaMEM executable
lamem = '"/home/chris/software/LaMEM/bin/opt/LaMEM"';
% LaMEM input file
input = '"Subduction2D.dat"';


if not(isfolder('Matrices'))
    mkdir('Matrices');
else
    rmdir('Matrices','s');
    mkdir('Matrices');
end

nel_x = 32;
nel_y = 2;
nel_z = 32;

% %% extract decomposition matrices
preMat = extract_preMat(lamem, input, nel_x,nel_y,nel_z);
preMat(:,:,end) = [];  % last matrix was only built for scaling --> dump it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_p     =  nel_x*nel_y*nel_z;
n_vis   = ((nel_x+1)*nel_y*nel_z)+(nel_x*(nel_y+1)*nel_z)+(nel_x*nel_y*(nel_z+1));
n_tot = n_vis + n_p;

N       = n_p+((nel_x+1)*(nel_y+1)*nel_z)+((nel_x+1)*nel_y*(nel_z+1))+(nel_x*(nel_y+1)*(nel_z+1));


%% dumped matrix before solving the system

%setup2D(par1(k1),par2(k2),par3(1));
% run simulation
system([lamem,' -ParamFile Subduction2D.dat -only_matrix']);

A_bs   =  PetscBinaryRead('Matrices/Ass_A.bin');
M_bs   =  PetscBinaryRead('Matrices/Ass_M.bin');
rhs_bs =  -PetscBinaryRead('Matrices/rhs.bin');
eta    =  PetscBinaryRead('Matrices/eta.bin');
rho    =  PetscBinaryRead('Matrices/rho.bin');
zrows  =  PetscBinaryRead('Matrices/zrows.bin');
zrows  = zrows +1;

J_bs = A_bs-M_bs;
solo = J_bs\rhs_bs;
%A_bs   =  full(A_bs-M_bs);
A_bsv  = A_bs(1:n_vis,1:n_vis);



%% assembled decomposition matrix
A_ass = sparse(n_tot,n_tot);

    for i = 1:N  
        temp = eta(i)*sparse((preMat(:,:,i)));
        i
        % multiply by reduced basis
        A_ass = A_ass + temp;
    end
  A_ass = A_ass(1:n_vis,1:n_vis);
    
% add bc (one's)
   % boundary condition matrix
A_ass           = full(A_ass);   
A_bc            = zeros(n_vis,n_vis); 
diagonal        = zeros(1,n_vis);
diagonal(zrows) = 1;
D               = diag(diagonal);
A_bc            = A_bc + D;
A_ass           = A_ass +A_bc;


A_diff = abs(A_ass)-abs(A_bsv);
spy(A_diff);

indices = find(abs(A_diff)<(1e-14));
A_diff(indices) = 0;
figure(2)
spy(A_diff);
max_dump = max(max(abs(A_bsv)));
max_ass  = max(max(abs(A_ass)));
max_error = max(max(abs(A_diff)));
   
%% dumped matrix after the solving system

%setup2D(par1(k1),par2(k2),par3(1));
%run simulation
system([lamem,' -ParamFile Subduction2D.dat']);

A_as   =  sparse(PetscBinaryRead('Matrices/Ass_A.bin'));
M_as   =  sparse(PetscBinaryRead('Matrices/Ass_M.bin'));
rhs_as =  -PetscBinaryRead('Matrices/rhs.bin');
sol = PetscBinaryRead('Matrices/sol.bin');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate residual

% interpolate density

n_nz    = (nel_x*nel_y*(nel_z-1));
n_vy = (nel_x*(nel_y+1)*nel_z);
n_vx = ((nel_x+1)*nel_y*nel_z);
n_xy    = nel_x*nel_y;
rho_i = zeros(n_nz,1);

for i = 1:n_nz
    rho_i(i) = (rho(i) + rho(i+(n_xy)))/2;
end

%rhs_blank = rhs_bs./rho;
test = nonzeros(rhs_bs);
rhs_blank = test./rho_i;
rhs_b = zeros(n_tot,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% uniiiits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate solution for dumped matrices
J_bs = A_bs-M_bs;
sol_d = J_bs\rhs_bs;

J_as = J_bs;
J_as(1:n_vis,1:n_vis) = 0;
J_as(1:n_vis,1:n_vis) = A_ass;

rhs_as = rhs_blank.*rho_i;
rhs_b(((((nel_x+1)*nel_y*nel_z) + ((nel_y+1)*nel_x*nel_z) + (nel_y*nel_x))+1):(((nel_x+1)*nel_y*nel_z) + ((nel_y+1)*nel_x*nel_z) + (nel_y*nel_x))+length(rho_i)) = rhs_as;

sol_ass = J_as\rhs_b;

diff_d_as = abs(sol_d) - abs(sol_ass);
max(diff_d_as)

diff_d_s = abs(sol_d) - abs(sol);
max(diff_d_s)

diff_ass_s = abs(sol_ass) - abs(sol);
max(diff_ass_s)

%% results are identical!

solz = sol(((((nel_x+1)*nel_y*nel_z) + ((nel_y+1)*nel_x*nel_z) + (nel_y*nel_x))+1):(((nel_x+1)*nel_y*nel_z) + ((nel_y+1)*nel_x*nel_z) + (nel_y*nel_x))+length(rho_i));


%% conversion from output to cm/yr:
length = 1000;
time   = 1e11;
cm_yr  = 0.01/(3600*24*365.25);
fac    =length/time/cm_yr;

