% function that creates a reduced basis of parameter space with greedy algorithm

% input:
% nel_x --> number of elements
% par   --> parameter space
% tol   --> tolerance for greedy algorithm

% output:
% B     --> reduced basis

function [B, res_max, ETA, RHO] = Reduced_Basis (lamem, input, nel_x, nel_y, nel_z,par,tol,n)

tic
%initialization
err       = 100;
loc       = 1;    % choosing first paramter sample to be evaluated
B         = [];   % initialize basis B
ETA       = [];   % initialize basis eta
RHO       = [];   % initialize basis rho

it        = 1;
res_max   = []; % vector to store maximal errors

% total number of velocity nodes in the mesh
n_velx = (nel_x+1)*nel_y*nel_z;
n_vely = (nel_y+1)*nel_x*nel_z;
n_velz = (nel_z+1)*nel_x*nel_y;
n_vel  = n_velx + n_vely + n_velz; % total number of all velocity points


%% greedy algorithm
while err > tol
     
    disp(['basis no. ',num2str(it)]);
    it = it +1 ;

    %% create truth solution and add it to basis
 
    % create markers
    setup2D(par(loc,:));
    %run simulation
    [t1, t2] = system([lamem,' -ParamFile ', input]);
    % read solution 
    sol_lamem =  PetscBinaryRead('Matrices/sol.bin');
    eta       =  PetscBinaryRead('Matrices/eta.bin');
    rho       =  PetscBinaryRead('Matrices/rho.bin');
    
    B = [B sol_lamem];   % enrich reduced basis
    ETA = [ETA eta];     % enrich eta basis
    RHO = [RHO rho];     % enrich rho basis
    
    % orthogonalization every n steps 
     if mod(it,n) == 0
     B = orth(B);
     end

    res_vec  = [];     % clear residual vector

    
    %% evaluate argmax of parameter space
    textprogressbar('parameter loop: ');
    for k = 1:length(par)
       
        %disp(['parameter loop: ',num2str((k/length(par))*100),'%'])
        % progressbar for parameter loop
          kp = (k/length(par))*100;
          textprogressbar(kp);

            % create markers
             setup2D(par(k,:));
            % run simulation
            [t1,t2] = system([lamem,' -ParamFile ', input, ' -only_matrix']);
            
            A   =  sparse(PetscBinaryRead('Matrices/Ass_A.bin'));
            M   =  sparse(PetscBinaryRead('Matrices/Ass_M.bin'));
            rhs =  -PetscBinaryRead('Matrices/rhs.bin');
                    
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
            %% K only with viscosity matrixes!!!!
            f = B.' * rhs;
                         
            alpha = K\f;
            Sol_RB = B * alpha;

            % access error --> global residual is used as an error estimator
            res = rhs(1:n_vel) - (VV * Sol_RB(1:n_vel)) - (PV * Sol_RB(n_vel+1:end));
            % maximum absolute value of res (maybe other norm??)
            maxres = max(abs(res));
            
            %disp('**********************');
            %disp(['res = ', num2str(maxres)]);
            %disp('**********************');
            res_vec = [res_vec maxres];   % vector of residual norm
                        
    end
    textprogressbar('');   
    [err, loc] = max(res_vec); % store max. error and location of max. error
    
    disp('*********************************************');
    disp(['max. error ',num2str(err), ' for parameter set no. ',num2str(loc)]);
    disp('*********************************************');
    res_max = [res_max err];


end

time = toc;
disp('==================================================================================');
disp('==============================   Greedy algorithm   ==============================');
disp('==================================================================================');
disp(['greedy algorithm is finished!!']);
disp(['duration: ',num2str(time),' s']);
disp(['total number of basis funtions: ',num2str(it)]);
disp('==================================================================================');
disp('==============================   Greedy algorithm   ==============================');
disp('==================================================================================');
end
