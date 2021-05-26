% precompute matrixes for reduced basis application
function [PMats] = precompute_mat (B, nel_x,nel_y,nel_z, mode, U)

%% initialization
% dimension of reduced basis
m       = length(B(1,:));
% number of pressure points / number of nodes
n_p     =  nel_x*nel_y*nel_z;   % number of pressure nodes
% dimension of velocity matrices
n_vis   = ((nel_x+1)*nel_y*nel_z)+(nel_x*(nel_y+1)*nel_z)+(nel_x*nel_y*(nel_z+1));
% number of dumped decomposition matrices
N       = n_p+((nel_x+1)*(nel_y+1)*nel_z)+((nel_x+1)*nel_y*(nel_z+1))+(nel_x*(nel_y+1)*(nel_z+1));
% context to store matrices

%% collect all decompositionj matrices (usually many --> N!!)
system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(1),' -eta[1] ', num2str(10),' -dump_decomposition_matrices']);

% assemble matrix
zrows   =  sparse(PetscBinaryRead('Matrices/zrows.bin'));
zrows   =  zrows+1; % Petsc starting indexing at 0, matlab at 1!

A   =  PetscBinaryRead('Matrices/Ass_A.bin');
divV    = sparse(A(end-n_p+1:end,1:end-n_p));
divP    = sparse(A(1:end-n_p,end-n_p+1:end));
B_t     =  B.';


%% without DEIM
if mode == 0
    PMats   = zeros(m,m,N);
    for i = 1:N
        
        A_ass = sparse(n_vis+n_p,n_vis+n_p);% extract pressure and velocity divergence matrices from assembled A matrix
        
        % load decomposition matrices and multiply by corresponding eta value
        Vis_fix = sparse(PetscBinaryRead(['Matrices/Vis',num2str(i),'.bin']));
        
        % boundary conditions
        Vis_fix(zrows,:) = 0;
        diagonal         = sparse(1,n_vis);
        diagonal(zrows)  = 1/N;
        D                = diag(diagonal);
        Vis_fix          = Vis_fix + D;
        
        A_ass(1:n_vis,1:n_vis) = Vis_fix;
        A_ass(end-n_p+1:end,1:end-n_p) = divV;
        A_ass(1:end-n_p,end-n_p+1:end) = divP;

        % multiply by reduced basis
        PMats(:,:,i) = B_t * A_ass * B;
        i
    end
    
%% with DEIM
elseif mode == 1
        dim_DEIM = length(U(1,:));
        PMats    = zeros(m,m,dim_DEIM);
        
    for j = 1:dim_DEIM
                
        Mi = sparse(m,m);
                
        for i = 1:N
                    
            A_ass = sparse(n_vis+n_p,n_vis+n_p);% extract pressure and velocity divergence matrices from assembled A matrix

            % load decomposition matrices and multiply by corresponding eta value
            Vis_fix = sparse(PetscBinaryRead(['Matrices/Vis',num2str(i),'.bin']));

            % boundary conditions
            Vis_fix(zrows,:) = 0;
            diagonal         = sparse(1,n_vis);
            diagonal(zrows)  = 1/N;
            D                = diag(diagonal);
            Vis_fix          = Vis_fix + D;

            A_ass(1:n_vis,1:n_vis) = Vis_fix;
            A_ass(end-n_p+1:end,1:end-n_p) = divV;
            A_ass(1:end-n_p,end-n_p+1:end) = divP;

            Mi = Mi + U(i,j)*(B_t * A_ass * B);


            % multiply by reduced basis
            i
            j
         end
            PMats(:,:,j) = Mi; 
    end
        
end
end