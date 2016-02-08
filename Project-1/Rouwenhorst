function P = Rouwenhorst(p, q, N)
%Rouwenhorst method: an alternative to Tauchen's method;
%P is a NxN transition matrix
%inputs:
%p, q: parameters, for AR(1) process, we can set p = q = (1+rho)/2
%N: dimension of matrix P
P = [p 1-p; 1-q, q];

for i = 3:N;
    P = p*[P zeros(i-1,1); zeros(1, i-1) 0]+(1-p)*[zeros(i-1,1) P; 0 zeros(1, i-1)] +(1-q)*[zeros(1, i-1) 0; P zeros(i-1,1)]+q*[0 zeros(1, i-1);zeros(i-1,1) P];
    P = [P(1,:); zeros(i-1, i)] + [zeros(1, i); 0.5*P(2:i-1,:); zeros(1, i)] + [zeros(i-1, i); P(i,:)] ;
end;

end
