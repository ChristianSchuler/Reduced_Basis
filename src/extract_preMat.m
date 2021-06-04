% extracts all decomposition Matrices and dumps a 3D sparse element
function [preMat] = extract_preMat (lamem, input, nel_x,nel_y,nel_z)

tic
% number of pressure points / number of nodes
n_p     =  nel_x*nel_y*nel_z;   % number of pressure nodes
% number of dumped decomposition matrices
N       = n_p+((nel_x+1)*(nel_y+1)*nel_z)+((nel_x+1)*nel_y*(nel_z+1))+(nel_x*(nel_y+1)*(nel_z+1));


[t1, t2] = system([lamem,' -ParamFile ', input, ' -dump_decomposition_matrices']);

% assemble matrix
zrows   =  sparse(PetscBinaryRead('Matrices/zrows.bin'));
zrows   =  zrows+1; % Petsc starting indexing at 0, matlab at 1!

points = [];
vals = [];
for i = 1:N
    

    % load decomposition matrices and multiply by corresponding eta value
    preMat = sparse(PetscBinaryRead(['Matrices/Vis',num2str(i),'.bin']));
    [k,j,s] = find(preMat);
    nump = length(k);
    i
    k
    for m = 1:nump
          
        if (ismember(k(m),zrows) == 0)  
        
         points = [points; k(m) j(m) i];
         vals   = [vals s(m)];
         
        else
            
         points = [points; k(m) j(m) i];
         vals   = [vals 0];
            
        end
         
    end      

end

preMat = ndSparse.build(points ,vals);

disp('extract dumped matrices and put it in sparse 3D structure:'); 
toc
end