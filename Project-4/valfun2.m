function val=valfun2(k)
global v0 beta delta alpha kmat k0 sigma kgrid P logzgrid

klo=max(sum(k>kmat),1); % identify the gridpoint that falls just below . .
% . . the choice for k
khi=klo+1;

% do the interpolation
gg = v0(klo,:) + (k-kmat(klo))*(v0(khi,:) - v0(klo,:))/(kmat(khi) - kmat(klo));
%gg = interp1(kmat,v0,k,'linear');
%exp(exp(logzgrid(2)))=2.7183
c = 2.7183*k0^alpha - k + (1-delta)*k0; % consumption
if c<=0
val = -8888888888888888-800*abs(c);
else
val = (1/(1-sigma))*(c^(1-sigma)-1) + beta*(gg*P(2,:)');
end
val = -val; % make it negative since we're maximizing and code is to minimize.
