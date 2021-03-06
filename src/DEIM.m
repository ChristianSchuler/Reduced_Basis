% DEIM - Discrete Emperical Interpolation Method

%input
% B - reduced basis of problem, obtained e.g by greedy algorithm
function [U,P,p] = DEIM (B)

% max of first basis vector
[val, loc] = max(abs(B(:,1)));

%initialize U,P and p

r      = length(B(1,:));   % rank
U      = [B(:,1)];
P      = zeros(length(B(:,1)),1);
P(loc) = 1;
p      = loc;

    for i = 2:r 
        % solve for c
        c           = (P.' * U) \ (P.' * B(:,i));
        %residual
        res         = B(:,i) - U*(c);
        [val, loc]  = max(abs(res));
        % update U,P and p
        U           = [U, B(:,i)];
        u_vec       = zeros(length(B(:,1)),1);
        u_vec(loc)  = 1;
        P           = [P, u_vec];
        p           = [p; loc]; 
    end
end