var K C L V U w q Y M eta mu b;
varexo e;
parameters theta alpha beta chi gamma delta_k rho sigma T delta_l G Phi A;
theta =0.7;
alpha= 0.34;
beta = 0.99;
chi = 1;
gamma =0.5;
delta_k =0.025;
rho= 0.80;
sigma = 0.1;
T = 1;
delta_l =  0.25;
G =0;
Phi =0.5;
A=1;

model;
L = (1 - delta_l)*L(-1)+M;
M = exp(b)*(U^theta)*V^(1-theta);
C+(1/T)*(K-(1-delta_k)*K(-1))+G = Y;
Y = A*(K^alpha)*L^(1-alpha);
1/T = beta*(C/C(+1))*A*alpha*(L(+1)/K(+1))^(1-alpha) + 1/T*(1-delta_k);
w -C*chi*(L+V+U)^gamma = mu - beta*(C/C(+1))*mu(+1)*(1-delta_l(+1));
A*(1-alpha)*(K/L)^(alpha)-w = eta -beta*(C/C(+1))*eta(+1)*(1-delta_l(+1));
-C*chi*(L+V+U)^gamma+ q =0;
eta*M/V - q = 0;
C*chi*(L+V+U)^gamma = mu*M/U;
eta = (1-Phi)*(eta+mu);
b = rho*b(-1)+e;
end;


initval;
K =1.0954;
C =0.2317;
L = 0.1233; 
V = 0.0308;
U = 0.0308;
w = 3;
q = 0.0996;
eta = 0.0996;
mu = 0.0996;
b=0;
end;


steady;

check;

shocks;
var e = sigma^2;
end;

stoch_simul(periods=2100);
