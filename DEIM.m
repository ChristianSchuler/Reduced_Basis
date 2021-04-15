% Creating interpolation of stiffness matrix with 'discrete imperical
% interploation method (DEIM)'

%input
% B - reduced basis of problem, obtained e.g by greedy algorithm

% output
% p - interpolation indices
function [P,p] = DEIM (B)

% max of first basis vector
[val, loc] = max(abs(B(:,1)));

%initialize U,P and p
U      = [B(:,1)];
P      = zeros(length(B(:,1)),1);
P(loc) = 1;
p      = loc;

    for i = 2:(length(B(1,:)))
        
        % solve for c
        c           = (P.' * B(:,i)) \ (P.' * U);
        
        %residual
        res         = B(:,i) - U*(c.');
        [val, loc]  = max(abs(res));
        
        % update U,P and p
        U           = [U, B(:,i)];
        u_vec       = zeros(length(B(:,1)),1);
        u_vec(loc)  = 1;
        P           = [P, u_vec];
        p           = [p; loc];
        
    end

end