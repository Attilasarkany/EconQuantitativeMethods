%%Value Function Iteration
%V(k,z) = u(c) + beta*sum(Pzz'*V(k',z'))
%c + k' = (e^z)*(k^alpha)+(1-delta)*k
%log(z') = rho*log(z) + sigma*epsilon
%epsilon~N(0,1)

clear all;clc; close all;
%% Parameter values
alpha = 0.66; % %parameter for the production function
beta = 0.8; %discout factor
delta = 0.025; %depreciation rate
rho = 0.9; %paramter for the AR process
sigma = 2; %parameter for the utility function
N = 1000; %number of grid points

%% Solve for steady state
kstar = ((1/beta - 1 + delta)/(exp(1)*alpha))^(1/(alpha-1));

%% Discretize the state space
kgrid = linspace(0.75*kstar, 1.25*kstar, N)';

logzgrid = [-0.05 0 0.05]';
nstate = length(logzgrid);

%% transition matrix
%Rouwenhorst method: an alternative to Tauchen's method;
%P is a NxN transition matrix
%inputs:
%p, q: parameters, for AR(1) process, we can set p = q = (1+rho)/2
%N: dimension of matrix P
p = (1+rho)/2;
q = (1+rho)/2;
P = Rouwenhorst(p, q, nstate); % P is the transition matrix for log(z)


%% Set up return matrix
U = zeros(N,N);
C = zeros(N,N);

for j = 1:nstate
for i = 1:N,
    k = kgrid(i);
    for h = 1:N,
        kprime = kgrid(h);
        c = exp(exp(logzgrid(j)))*k^alpha + (1-delta)*k - kprime;
        %k, kprime, logz
        C(i,h,j) = max(0,c); %% Imposes consumption non-negative
        U(i,h,j) = (1/(1-sigma))*(max(0,c)^(1-sigma)-1);
    end;
end;
end;

%% Begin value function iteration
iter = 0;
error = 100;
V0 = zeros(N,nstate);

while error>10^(-3) && iter<1000        
   for j = 1:nstate
       %given k, logz, find the best kprime
    TV1 = max((U(:,:,j)+beta*ones(N,1)*P(j,:)*V0')'); %% Bellman operation   
    TV1 = TV1';
    TV0(:,j) = TV1;
   end 
error = norm(V0-TV0);
V0 = TV0; %% If not, the loop updates; otherwise, it stops
iter = iter+1
end
V = TV0;


%%
figure 
plot(kgrid,V(:,1))
hold on
plot(kgrid,V(:,2))
plot(kgrid,V(:,3))
hold off
legend('V(log(z)=-0.05)', 'V(log(z)= 0 )', 'V(log(z)=0.05)', 2, 'Location','northwest')
title('Value function')
xlabel('Capital stock')
ylabel('Lifetime utility')


%% Compute policy function 

gs = zeros(N,3);
cs = zeros(N,3);
for j = 1:nstate
   [vs,is] = max((U(:,:,j)+beta*ones(N,1)*P(j,:)*V')');
   for i = 1:N,
   gs(i,j) = kgrid(is(i));
   cs(i,j) = (1/(1-sigma))*(exp(exp(logzgrid(j)))*(kgrid(i)^alpha+(1-delta)*kgrid(i)-kgrid(is(i)))^(1-sigma)-1);
   end
end

figure 
plot(kgrid,kgrid,'--k')
hold on
plot(kgrid,gs(:,1))
plot(kgrid,gs(:,2))
plot(kgrid,gs(:,3))
hold off
legend('45-degree','k_{t+1}=g(k_t, log(z)_t=-0.05)', 'k_{t+1}=g(k_t, log(z)_t=0)', 'k_{t+1}=g(k_t, log(z)_t=0.05)', 'Location','northwest')
title('Policy function')
xlabel('k_{t}')
ylabel('k_{t+1}')
