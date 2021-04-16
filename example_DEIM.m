clear all,clc

% nonlinear test function
test_func = @(x,p) ((1-x).*cos(3*pi*p*(x+1))).*exp(-(1+x)*p);
x = linspace(-1,1,100);

%parameter space
p_pick = @(st,en) st + rand*(en-st);

% create snapshots
n = 50;
Sn = [];
for i = 1:n
    a = test_func(x,p_pick(1,pi));
    Sn = [Sn, a.']; % snapshot matrix
end

% rank
r = 10;
% POD/SVD
[U,S,V] = svd(Sn,0);

% choose basis functions with rank r
U = U(:,1:r);

[U,P,p] = DEIM(U);  % applying DEIM function
  
% calculate DEIM approximation of solution of specific parameter
sol_int = (U*inv((P.'*U)) * P.')* test_func(x,1.5).';


plot(x,test_func(x,1.5),x,sol_int,'--','LineWidth',2);
title(['analytical solution vs. DEIM solution with POD rank ', num2str(r)]);
legend('analytical','DEIM');

