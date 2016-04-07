%V(k,z) = u(c) + beta*sum(Pzz'*V(k',z'))
%c + k' = (e^z)*(k^alpha)+(1-delta)*k
%log(z') = rho*log(z) + sigma*epsilon
%epsilon~N(0,1)
clear all;clc; close all;
%%
alpha = 0.66;
beta = 0.8;
delta = 
rho = 
sigma = 

%%
f.prod = @(k, z)(e^z)*(k^alpha);
u.log = @(c)(log(c));

initials = [13.079 0.3443];  %[ k z];
kmin = 0.75*initials(1); kmax = 1.25*initials(1);
zmin = 0.75*initials(2); zmax = 1.25*initials(2);
N = 50; 
kgrid = linspace(kmin, kmax, N);
zgrid = linspace(zmin, zmax, N);

[K, Z] = meshgrid(kgrid, zgrid);
