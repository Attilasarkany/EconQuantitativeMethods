clear all;clc; close all;
%%Reference: http://www3.nd.edu/~esims1/val_fun_iter.pdf
%Eric Sims, University of Notre Dame, Spring 2011
%Continuous time value function iteration using linear splines

%%
tic

global v0 beta delta alpha kmat k0 sigma kgrid P logzgrid

%% Parameter values
alpha = 0.66; % %parameter for the production function
beta = 0.8; %discout factor
delta = 0.025; %depreciation rate
rho = 0.9; %paramter for the AR process
sigma = 2; %parameter for the utility function
kgrid = 1000; %number of grid points

tol = 0.001;
maxits = 300;
dif = tol+1000;
its = 0;
%% Solve for steady state
kstar = ((1/beta - 1 + delta)/(exp(1)*alpha))^(1/(alpha-1));

%% Discretize the state space
%kgrid = linspace(0.75*kstar, 1.25*kstar, N)';

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

%%
kmin = 0.75*kstar; % minimum
kmax = 1.25*kstar; % maximum
grid = (kmax-kmin)/kgrid; % grid

kmat = kmin:grid:kmax;

kmat = kmat';

[N,n] = size(kmat);

v0 = zeros(N,3);
%%
v1 = zeros(N,3);
k11 = zeros(N,3);
while dif>tol && its < maxits
    for j = 1
        for i = 1:N
            k0 = kmat(i,1);
            a0 = logzgrid(j,1);
            k1 = fminbnd(@valfun1,kmin,kmax);
            v1(i,j) = -valfun1(k1);
            k11(i,j) = k1;
        end
    end
    for j = 2
        for i = 1:N
            k0 = kmat(i,1);
            a0 = logzgrid(j,1);
            k1 = fminbnd(@valfun2,kmin,kmax);
            v1(i,j) = -valfun2(k1);
            k11(i,j) = k1;
        end
    end
    for j = 3
        for i = 1:N
            k0 = kmat(i,1);
            a0 = logzgrid(j,1);
            k1 = fminbnd(@valfun3,kmin,kmax);
            v1(i,j) = -valfun3(k1);
            k11(i,j) = k1;
        end
    end
    dif = norm(v1 - v0)
    v0 = v1;
    its = its+1
end
toc
%%
figure 
plot(kmat,v1(:,1))
hold on
plot(kmat,v1(:,2))
plot(kmat,v1(:,3))
hold off
legend('V(log(z)=-0.05)', 'V(log(z)= 0 )', 'V(log(z)=0.05)', 2, 'Location','northwest')
title('Value function')
xlabel('Capital stock')
ylabel('Lifetime utility')

