clear all,clc

% see chaturantabut 2010 sec. 3.3.1 for reference

% nonlinear test function
test_func = @(x,p) ((1-x).*cos(3*pi*p*(x+1))).*exp(-(1+x)*p);
x = linspace(-1,1,100);

%parameter space
%p_pick = @(st,en) st + rand*(en-st);

% create snapshots
n = 51;
Sn = [];
par = linspace(1,pi,51);
for i = 1:n
    %a = test_func(x,p_pick(1,pi));
    a = test_func(x,par(i));
    Sn = [Sn, a.']; % snapshot matrix
end

% rank
r = 10;
% POD/SVD
[U,S,V] = svd(Sn,0);


% choose basis functions with rank r
U = U(:,1:r);

[U,P] = DEIM(U);      % applying DEIM function
[row, col] = find(P);
ip = row;             % interpolation points

%plot first six POD bases
for i= 1:6
figure(1)
plot(x,U(:,i));
hold on
end
legend('POD basis 1','POD basis 2','POD basis 3','POD basis 4','POD basis 5','POD basis 6');

ni = 3.1;
% calculate DEIM approximation of solution of specific parameter
i_m    = inv(P.'*U) * P.';         % interpolation matrix
ct      = i_m * test_func(x,ni).';
sol_int = U*ct;   % test_func(x,ni) (=ungefaehr) U*ct

figure(2)
plot(x,test_func(x,ni),x,sol_int,'--','LineWidth',2);
hold on;
scatter(x(ip),test_func(x(ip),ni),'filled','red');
title(['analytical solution vs. DEIM solution with POD rank ', num2str(r)]);
legend('analytical','DEIM','DEIM interpolation points');

