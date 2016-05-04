function val=valfun1(k)
%%Reference: http://www3.nd.edu/~esims1/val_fun_iter.pdf
%Eric Sims, University of Notre Dame, Spring 2011
%Continuous time value function iteration using linear splines

%%
global v0 beta delta alpha kmat k0 sigma  P 
kmin=max(sum(k>kmat),1); 
kmax=kmin+1;

% do the interpolation
gg = v0(kmin,:) + (k-kmat(kmin))*(v0(kmax,:) - v0(kmin,:))/(kmat(kmax) - kmat(kmin));

% exp(exp(logzgrid(1,1)))
c =2.5889*k0^alpha - k + (1-delta)*k0; % consumption
if c<=0
val = -8888888888888888-800*abs(c);
else
val = (1/(1-sigma))*(c^(1-sigma)-1) + beta*(gg*P(1,:)');
end
val = -val; % make it negative since we're maximizing and code is to minimize.