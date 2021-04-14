% input:
% nel_x --> number of elements
% par   --> parameter space
% tol   --> tolerance for greedy algorithm

% output:
% B     --> reduced basis

function [B,res,adj_vec, res_max] = create_RB (nel_x, nel_y, nel_z,par1,par2,tol,method)

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
[temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(loc1)),' -eta[1] ', num2str(par2(loc2)), ' -eta[2] ', num2str(par2(loc2)),' -eta[3] ', num2str(par2(loc2)),' -eta[4] ', num2str(par2(loc2)),' -eta[5] ', num2str(par2(loc2)),' -eta[6] ', num2str(par2(loc2)), ' -eta[7] ', num2str(par2(loc2)),' -eta[8] ', num2str(par2(loc2)),' -eta[9] ', num2str(par2(loc2)),' -eta[10] ', num2str(par2(loc2))]);
%[temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ./FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(loc1)),' -rho[1] ', num2str(par2(loc2)), ' -rho[2] ', num2str(par2(loc2)),' -rho[3] ', num2str(par2(loc2)),' -rho[4] ', num2str(par2(loc2)),' -rho[5] ', num2str(par2(loc2)),' -rho[6] ', num2str(par2(loc2)), ' -rho[7] ', num2str(par2(loc2)),' -rho[8] ', num2str(par2(loc2)),' -rho[9] ', num2str(par2(loc2)),' -rho[10] ', num2str(par2(loc2))]);

 % read data 
A   =  PetscBinaryRead('Mono_A.bin');
M   =  PetscBinaryRead('Mono_M.bin');
rhs =  PetscBinaryRead('r.1.bin');
J = A;% - M;

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
            [temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ../FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(k1)),' -eta[1] ', num2str(par2(k2)), ' -eta[2] ', num2str(par2(k2)),' -eta[3] ', num2str(par2(k2)),' -eta[4] ', num2str(par2(k2)),' -eta[5] ', num2str(par2(k2)),' -eta[6] ', num2str(par2(k2)), ' -eta[7] ', num2str(par2(k2)),' -eta[8] ', num2str(par2(k2)),' -eta[9] ', num2str(par2(k2)),' -eta[10] ', num2str(par2(k2))]);
            %[temp1, temp2] = system(['/home/chris/software/LaMEM/bin/opt/LaMEM -ParamFile ./FallingBlock_mono_PenaltyDirect.dat -eta[0] ', num2str(par1(loc1)),' -rho[1] ', num2str(par2(loc2)), ' -rho[2] ', num2str(par2(loc2)),' -rho[3] ', num2str(par2(loc2)),' -rho[4] ', num2str(par2(loc2)),' -rho[5] ', num2str(par2(loc2)),' -rho[6] ', num2str(par2(loc2)), ' -rho[7] ', num2str(par2(loc2)),' -rho[8] ', num2str(par2(loc2)),' -rho[9] ', num2str(par2(loc2)),' -rho[10] ', num2str(par2(loc2))]);

            % read data 
            A   =  PetscBinaryRead('Mono_A.bin');
            M   =  PetscBinaryRead('Mono_M.bin');
            rhs =  PetscBinaryRead('r.1.bin');
            J = A;% - M;

            if method == 1

                %==========================================================
                % create RB with residual over whole domain
                %==========================================================

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

            elseif method == 2

                %=============================================================
                % create RB with the help of the adjoint
                %=============================================================
               
                % Quantitiy of interest 
                % rhs vector
                rhs_adj = zeros(length(Sol_T),1);
                %rhs_adj(4352*2:4352*3) = 1;
                %rhs_adj(4352:4352*2) = 1;
                rhs_adj(1:4352) = 1;
                %rhs_adj(4352*3:end) = 1;


                % solve adjoint system
                adj_vec = (J)\rhs_adj;   % J is symmetric --> J^T = J
                
                % create RB solution
                K = B.' * (J) * B;
                f = B.' * rhs;
                alpha = K\f;
                Sol_RB = B * alpha;

                % access error --> global residual is used as an error estimator
                res = adj_vec.' * (rhs - (J*Sol_RB));
                disp(['**********************']);
                disp(['res = ', num2str(abs(res))]);
                disp(['**********************']);
    %            temp_res = (rhs - (J*Sol_RB));
    %            res = dot(adj_vec,temp_res)
                res_vec = [res_vec res];   % irgendeine coole norm über das residual

            end 
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
