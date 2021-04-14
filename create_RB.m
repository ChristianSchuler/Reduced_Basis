% input:
% nel_x --> number of elements
% par   --> parameter space
% tol   --> tolerance for greedy algorithm

% output:
% B     --> reduced basis

function [B,res,adj_vec, res_max] = create_RB (nel_x, nel_y, nel_z,par1,par2,tol)

%initialization
err      = 100;
loc1      = 1;    % choosing first paramter of par as first guess
loc2      = 1;    % choosing first paramter of par as first guess
B        = [];   % initialize basis B
it       = 1;
res_max  = []; % vector to store maximal errors

%% greedy algorithm
while err > tol
    
disp(['basis no. ',num2str(it)]);
it = it +1 ;

%% create truth solution and add it to basis
% create Jacobian and rhs vector solution
[temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(loc1)),' -eta[1] ', num2str(par1(loc1))]);
%[temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(loc1)),' -rho[1] ', num2str(par2(loc2))]);

 % read data 
A   =  PetscBinaryRead('Mono_A.bin');
M   =  PetscBinaryRead('Mono_M.bin');
rhs =  PetscBinaryRead('r.1.bin');

% creating truth solution
[Sol_T,Sol_Vel,Sol_P,VV,VP,PV,PP] = solve_stokes(A,M,rhs,nel_x, nel_y, nel_z);

B = [B Sol_T];              % add solution to reduced basis
B = orth(B);

res_vec  = []; % empty residual vector
 
%% evaluate argmax of parameter space
it2 = 0;
it3 = 0;
    for k1 = 1:length(par1)
        
        %disp(['parameter loop: ',num2str(it2+it3/(length(par))*100),'%']);
        %it2 = it2 +1 ;
        
        for k2 = 1:length(par2)
            
            it3 = it3 +1 ;
            disp(['parameter loop: ',num2str(((it3)/((length(par1)*(length(par2)))))*100),'%']);
            
            
            % create Jacobian and rhs vector solution
            [temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(loc1)),' -eta[1] ', num2str(par1(loc1))]);
            %[temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(loc1)),' -rho[1] ', num2str(par2(loc2))]);

            % read data 
            A   =  PetscBinaryRead('Mono_A.bin');
            M   =  PetscBinaryRead('Mono_M.bin');
            rhs =  PetscBinaryRead('r.1.bin');
            J = A;% - M;

                %=============================================================
                % create RB with residual over whole domain
                %=============================================================

                % solve RB
                K = B.' * (J) * B;
                f = B.' * rhs;
                alpha = K\f;
                Sol_RB = B * alpha;

                % access error --> global residual is used as an error estimator
                res = rhs - (J * Sol_RB);            % calculate residual
                disp(['**********************']);
                disp(['res = ', num2str(max(abs(res)))]);
                disp(['**********************']);
                res_vec = [res_vec max(abs(res))];   % irgendeine coole norm über das residual

                adj_vec = 0;
                
        end

    end
    
[err, loc] = max(abs(res_vec)); % store max. error and location of max. error
disp(['==========================================================']);
disp(['max. error ',num2str(err), ' for parameter no. ',num2str(loc)]);
disp(['==========================================================']);
res_max = [res_max err];

% parameter estimation for next truth solution
loc;
loc2 = mod(loc,length(par2));
if loc2 == 0
   loc2 = length(par2);
end
loc1 = (loc -loc2) / length(par2) + 1;


end

disp(['==============================================']);
disp(['finished!! total number of basis funtions ',num2str(it)]);
disp(['==============================================']);
