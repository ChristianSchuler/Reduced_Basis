% function that creates a reduced basis of parameter space with greedy algorithm

% input:
% nel_x --> number of elements
% par   --> parameter space
% tol   --> tolerance for greedy algorithm

% output:
% B     --> reduced basis

function [B, res_max, ETA, RHO] = Reduced_Basis (lamem, input, nel_x, nel_y, nel_z,g,par1,par2,tol)

%initialization
err       = 100;
loc1      = 1;    % choosing first paramter of par as first guess
loc2      = 1;    % choosing first paramter of par as first guess
B         = [];   % initialize basis B
ETA       = [];   % initialize basis eta
RHO       = [];   % initialize basis rho
it        = 0;
res_max   = []; % vector to store maximal errors


%% variables for calculating rhs
n_nz    = (nel_x*nel_y*(nel_z-1));
n_xy    = nel_x*nel_y;

 % total number of velocity nodes in the mesh
n_velx = (nel_x+1)*nel_y*nel_z;
n_vely = (nel_y+1)*nel_x*nel_z;
n_velz = (nel_z+1)*nel_x*nel_y;
n_vel  = n_velx + n_vely + n_velz; % total number of all velocity points

% pressure nodes
n_p = nel_x * nel_y * nel_z;

%total nodes
N = n_vel+n_p;

% rhs wo rho
g_fac = g/2;

inputG      = input;
inputG(1)   = [];
inputG(end) = [];


%% greedy algorithm
while err > tol
     
    disp(['basis no. ',num2str(it)]);
    it = it +1 ;

    %% create truth solution and add it to basis
    % create Jacobian and rhs vector solution
    %copyfile(['../',inputG],'geometry.dat');
%     copyfile(inputG,'geometry.dat');
%     radius 		= 	par1(loc1);	
%     str = ['<SphereStart>' newline 'phase  = 1' newline 'center = 0.5 0.5 0.5' newline 'radius = ' num2str(radius) newline '<SphereEnd>'];
%     fid=fopen('geometry.dat','a+');
%     fprintf(fid, str);
%     fclose(fid);

%      [t1, t2] = system([lamem,' -ParamFile geometry.dat']);
    %[t1, t2] = system([lamem,' -ParamFile ../', input, ' -eta[0] ', num2str(par1(loc1)),' -eta[1] ', num2str(par2(loc2))]);
    
    
    %create partitioning file
    %system(['mpiexec -n 4 ' , lamem,' -ParamFile Subduction3D.dat -mode save_grid']);
    % create markers
    setup3D(par1(loc1),par2(loc2));
    % run simulation
    [t1,t2] = system([lamem,' -ParamFile Subduction3D.dat']);
    
    % read data 
    
    %eta       =  PetscBinaryRead('Matrices/eta.bin');
    %rho       =  PetscBinaryRead('Matrices/rho.bin');
    sol_lamem =  PetscBinaryRead('Matrices/sol.bin');
    
    %ETA = [ETA eta];   % enrich eta basis
    %RHO = [RHO rho];   % enrich rho basis
    
    B = [B sol_lamem];     % enrich reduced basis
    
%     if mod(it,2) == 0
%     B = orth(B);
%     end

    res_vec  = [];     % clear residual vector

    
    %% evaluate argmax of parameter space
    it2 = 0;
    for k1 = 1:length(par1)


        for k2 = 1:length(par2)

            it2 = it2 +1 ;
            disp(['parameter loop: ',num2str(((it2)/((length(par1)*(length(par2)))))*100),'%']);
 
            % create 'Matrices' folder or first delete it if it already exists
            if not(isfolder('Matrices'))
            mkdir('Matrices');
            else
            rmdir('Matrices','s');
             mkdir('Matrices');
            end
                        
            % create Jacobian and rhs vector solution  
%             copyfile(inputG,'geometry.dat');
%             radius 		= 	par1(k1);	
%             str = ['<SphereStart>' newline 'phase  = 1' newline 'center = 0.5 0.5 0.5' newline 'radius = ' num2str(radius) newline '<SphereEnd>'];
%             fid=fopen('geometry.dat','a+');
%             fprintf(fid, str);
%             fclose(fid);
%             [t1, t2] = system([lamem,' -ParamFile geometry.dat -only_matrix']);
            %[t1, t2] = system([lamem,' -ParamFile ../', input, ' -eta[0] ', num2str(par1(k1)),' -eta[1] ', num2str(par2(k2)),' -only_matrix']);
            
            %create partitioning file
            %[t1,t2] = system(['mpiexec -n 4 ' , lamem,' -ParamFile Subduction3D.dat -mode save_grid']);
            % create markers
            setup3D(par1(k1),par2(k2));
            % run simulation
            [t1,t2] = system([lamem,' -ParamFile Subduction3D.dat -only_matrix']);
            
            A   =  sparse(PetscBinaryRead('Matrices/Ass_A.bin'));
            M   =  sparse(PetscBinaryRead('Matrices/Ass_M.bin'));
            rhs =  -PetscBinaryRead('Matrices/rhs.bin');
           
            
            if (it == 1)
            
                eta       =  PetscBinaryRead('Matrices/eta.bin');
                eta       =  eta*(1e21);

                rho =  PetscBinaryRead('Matrices/rho.bin');
                rho =  rho*(1e25);
                
                ETA = [ETA eta];   % enrich eta basis
                RHO = [RHO rho];   % enrich rho basis
            end
                  
            % extract Jacobian
            J   = A - M;
            
            % collect block matrixes
            VV = J(1:n_vel,1:n_vel);
            PV = J(1:n_vel,n_vel+1:end);

            %=============================================================
            % create RB with residual over whole domain
            %=============================================================

            % solve RB
            K = B.' * (J) * B;
            f = B.' * rhs;
            detK = det(K);
            
%             if detK < 1e-14
%                  maxres = 0;
%                  warning('matrix K is close to singular');
%               Bo = orth(B);
%               Ko = Bo.' * J * Bo;
%               fo = Bo.' * rhs;
%               alphao = Ko\fo;
%              Sol_RBo = Bo * alphao;

%             else
                
            alpha = K\f;
            Sol_RB = B * alpha;
     

             % access error --> global residual is used as an error estimator
            res = rhs(1:n_vel) - (VV * Sol_RB(1:n_vel)) - (PV * Sol_RB(n_vel+1:end));
            maxres = max(abs(res));
            %end
            
            disp(['**********************']);
            disp(['res = ', num2str(maxres)]);
            disp(['**********************']);
            res_vec = [res_vec maxres]; 
                      
        end
  
    end
    
    if max(res_vec) == 0
       break;
       error('reduced basis matrix is singular for every parameter');
    end
     
    [err, loc] = max(res_vec); % store max. error and location of max. error
    
    disp(['==========================================================']);
    disp(['max. error ',num2str(err), ' for parameter no. ',num2str(loc)]);
    disp(['==========================================================']);
    res_max = [res_max err];

    % parameter estimation for next truth solution
    loc2 = mod(loc,length(par2));
    if loc2 == 0
       loc2 = length(par2);
    end
    loc1 = (loc -loc2) / length(par2) + 1;


end

disp(['==================================================================================']);
disp(['reduced basis created with greedy algorithm!! total number of basis funtions ',num2str(it)]);
disp(['==================================================================================']);
end
