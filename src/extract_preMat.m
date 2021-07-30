% extracts all decomposition Matrices and dumps a 3D sparse element
function [preMat] = extract_preMat (lamem, input, nel_x,nel_y,nel_z)

tic
% number of pressure points / number of nodes
n_p     =  nel_x*nel_y*nel_z;   % number of pressure nodes
% dimension of velocity matrices
n_vis   = ((nel_x+1)*nel_y*nel_z)+(nel_x*(nel_y+1)*nel_z)+(nel_x*nel_y*(nel_z+1));
n_tot = n_vis + n_p;
% number of dumped decomposition matrices
N       = n_p+((nel_x+1)*(nel_y+1)*nel_z)+((nel_x+1)*nel_y*(nel_z+1))+(nel_x*(nel_y+1)*(nel_z+1));

disp('dumping decomposition matrices from LaMEM...')
[t1,t2] = system([lamem,' -ParamFile ', input, ' -dump_decomposition_matrices']);

% assemble matrix
zrows   =  sparse(PetscBinaryRead('Matrices/zrows.bin'));
zrows   =  zrows+1; % Petsc starting indexing at 0, matlab at 1!

points = [];
vals = [];

    textprogressbar('extracting decomposition matrices: ');
for i = 1:N 
   % load decomposition matrices and multiply by corresponding eta value
    preMat = sparse(PetscBinaryRead(['Matrices/Vis',num2str(i),'.bin']));
    [k,j,s] = find(preMat);
    nump = length(k);
   
    for m = 1:nump         
        if (ismember(k(m),zrows) == 0)          
         points = [points; k(m) j(m) i];
         vals   = [vals s(m)];     
        else     
         points = [points; k(m) j(m) i];
         vals   = [vals 0];   
         
        end      
    end  
    if (mod(i,100) == 0)
    ip = (i/N)*100;
    textprogressbar(ip);
    end
end
textprogressbar('');
    
%points = [points; n_tot n_tot N+1];
points = [points; n_vis n_vis N+1];
vals   = [vals 8];
preMat = ndSparse.build(points ,vals);

time = toc;
disp('=====================================================================');
disp('================  extracting decomposition matrices   ===============');
disp(['extract dumped matrices and put it in sparse 3D structure:']);
disp(['duration: ',num2str(time),' s']);
disp(['total number of dumped matrices: ',num2str(N)]);
disp('=====================================================================');

end