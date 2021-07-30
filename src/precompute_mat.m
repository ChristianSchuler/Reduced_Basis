% precompute matrixes for reduced basis application
function [PMats] = precompute_mat (preMat, B, nel_x,nel_y,nel_z, mode, U) %#codegen

%% initialization
% dimension of reduced basis
m       = length(B(1,:));
% number of pressure points / number of nodes
n_p     =  nel_x*nel_y*nel_z;   % number of pressure nodes
% dimension of velocity matrices
n_vis   = ((nel_x+1)*nel_y*nel_z)+(nel_x*(nel_y+1)*nel_z)+(nel_x*nel_y*(nel_z+1));
% number of dumped decomposition matrices
N       = n_p+((nel_x+1)*(nel_y+1)*nel_z)+((nel_x+1)*nel_y*(nel_z+1))+(nel_x*(nel_y+1)*(nel_z+1));
% one more matrix for (((((pressure velo divergence matrix/))))another one for
% boundary conditions
N       = N+1;

% assemble matrix
zrows   =  sparse(PetscBinaryRead('Matrices/zrows.bin'));
zrows   =  zrows+1; % Petsc starting indexing at 0, matlab at 1!

%A   =  PetscBinaryRead('Matrices/Ass_A.bin');
%divV    = sparse(A(end-n_p+1:end,1:end-n_p));
%divP    = sparse(A(1:end-n_p,end-n_p+1:end));
B   = B(1:n_vis,:);
B   = sparse(B);
B_t =  B.';


%% without DEIM
if mode == 0
    tic
    PMats   = zeros(m,m,N);
   
    
   disp('multiplying B with decomposition matrices: ');
    for i = 1:N-1  
        
       % multiply by reduced basis
        temp = (B_t * preMat(:,:,i)) * B;
        PMats(:,:,i) = temp;
        
    end
% pressure and velocity divergence submatrix
% A_ass                          = sparse(n_vis+n_p,n_vis+n_p);    
% A_ass(end-n_p+1:end,1:end-n_p) = divV;
% A_ass(1:end-n_p,end-n_p+1:end) = divP;
% PMats(:,:,end-1)               = B_t * A_ass * B;

% boundary condition matrix
%A_ass            = sparse(n_vis+n_p,n_vis+n_p); 
A_ass            = sparse(n_vis,n_vis); 
diagonal         = sparse(1,n_vis);
diagonal(zrows)  = 1;
D                = diag(diagonal);
A_ass            = A_ass + D;
PMats(:,:,end)   = B_t * A_ass * B;
 
time = toc;

disp('=====================================================================');
disp('============  multiply  matrices with B (without DEIM)  =============');
disp(['multiply viscosity matrices with reduced basis B']);
disp(['duration: ',num2str(time),' s']);
disp('=====================================================================');


%% with DEIM
elseif mode == 1
    tic
    U = sparse(U);
    dim_DEIM = length(U(1,:));
    PMats    = zeros(m,m,dim_DEIM+1);
        
    textprogressbar('multiplying B with decomposition matrices with DEIM: ');  
    for j = 1:dim_DEIM
                
        Mi = sparse(m,m);
         
        for i = 1:N-1    
                   
            temp = (B_t * preMat(:,:,i)) * B;

            % multiply with reduced basis 
            temp = U(i,j)*temp;
            Mi   = Mi + temp;  
            
        end
            PMats(:,:,j) = Mi; 
            
           jp = (j/dim_DEIM)*100;
           textprogressbar(jp);
    end
    textprogressbar('');
    
    disp('precompute matrix components with DEIM'); 
    time = toc;
    
% % pressure and velocity divergence submatrix
% A_ass                          = sparse(n_vis+n_p,n_vis+n_p);    
% A_ass(end-n_p+1:end,1:end-n_p) = divV;
% A_ass(1:end-n_p,end-n_p+1:end) = divP;
% PMats(:,:,end-1)               = B_t * A_ass * B;

% boundary condition matrix
A_ass            = sparse(n_vis,n_vis); 
diagonal         = sparse(1,n_vis);
diagonal(zrows)  = 1;
D                = diag(diagonal);
A_ass            = A_ass + D;
PMats(:,:,end)   = B_t * A_ass * B;

disp('=====================================================================');
disp('==============  multiply  matrices with B (with DEIM)  ==============');
disp(['multiply viscosity matrices with reduced basis B']);
disp(['duration: ',num2str(time),' s']);
disp('=====================================================================');


end


end