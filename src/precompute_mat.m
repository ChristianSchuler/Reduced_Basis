% precompute matrixes for reduced basis application
function [PMats] = precompute_mat (lamem, input, preMat, B, nel_x,nel_y,nel_z, mode, U) %#codegen

%% initialization
% dimension of reduced basis
m       = length(B(1,:));
% number of pressure points / number of nodes
n_p     =  nel_x*nel_y*nel_z;   % number of pressure nodes
% dimension of velocity matrices
n_vis   = ((nel_x+1)*nel_y*nel_z)+(nel_x*(nel_y+1)*nel_z)+(nel_x*nel_y*(nel_z+1));
% number of dumped decomposition matrices
N       = n_p+((nel_x+1)*(nel_y+1)*nel_z)+((nel_x+1)*nel_y*(nel_z+1))+(nel_x*(nel_y+1)*(nel_z+1));
% one more matrix for pressure velo divergence matrix/ another one for
% boundary conditions
N       = N+2;

%% collect all decompositionj matrices (usually many --> N!!)
[t1, t2] = system([lamem,' -ParamFile ', input, ' -dump_decomposition_matrices']);

% assemble matrix
zrows   =  sparse(PetscBinaryRead('Matrices/zrows.bin'));
zrows   =  zrows+1; % Petsc starting indexing at 0, matlab at 1!

A   =  PetscBinaryRead('Matrices/Ass_A.bin');
divV    = sparse(A(end-n_p+1:end,1:end-n_p));
divP    = sparse(A(1:end-n_p,end-n_p+1:end));
B_t     =  B.';
B_t     = sparse(B_t);
B       = sparse(B);




%% without DEIM
if mode == 0
    tic
    PMats   = zeros(m,m,N);
   
    for i = 1:N-2  
        
        temp = (B_t * preMat(:,:,i)) * B;
     
        % multiply by reduced basis
        PMats(:,:,i) = temp;
    end

% pressure and velocity divergence submatrix
A_ass                          = sparse(n_vis+n_p,n_vis+n_p);    
A_ass(end-n_p+1:end,1:end-n_p) = divV;
A_ass(1:end-n_p,end-n_p+1:end) = divP;
PMats(:,:,end-1)               = B_t * A_ass * B;

% boundary condition matrix
A_ass            = sparse(n_vis+n_p,n_vis+n_p); 
diagonal         = sparse(1,n_vis+n_p);
diagonal(zrows)  = 1;
D                = diag(diagonal);
A_ass            = A_ass + D;
PMats(:,:,end)   = B_t * A_ass * B;
    disp('precompute matrix components'); 
    toc


%% with DEIM
elseif mode == 1
    U = sparse(U);
    tic
    dim_DEIM = length(U(1,:));
    PMats    = zeros(m,m,dim_DEIM+2);
        
    for j = 1:dim_DEIM
                
        Mi = sparse(m,m);
                
        for i = 1:N-2    
                   
            temp = (B_t * preMat(:,:,i)) * B;

            % multiply with reduced basis 
            temp = U(i,j)*temp;
            Mi   = Mi + temp;   
         end
            PMats(:,:,j) = Mi; 
    end
    disp('precompute matrix components with DEIM'); 
    toc
    
% pressure and velocity divergence submatrix
A_ass                          = sparse(n_vis+n_p,n_vis+n_p);    
A_ass(end-n_p+1:end,1:end-n_p) = divV;
A_ass(1:end-n_p,end-n_p+1:end) = divP;
PMats(:,:,end-1)               = B_t * A_ass * B;

% boundary condition matrix
A_ass            = sparse(n_vis+n_p,n_vis+n_p); 
diagonal         = sparse(1,n_vis+n_p);
diagonal(zrows)  = 1;
D                = diag(diagonal);
A_ass            = A_ass + D;
PMats(:,:,end)   = B_t * A_ass * B;

end
end