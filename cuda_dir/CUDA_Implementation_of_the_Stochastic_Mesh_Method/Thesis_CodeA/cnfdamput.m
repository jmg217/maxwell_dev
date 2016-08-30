function put = cnfdamput(Smax, dS, T, dT, X, R, SIG)
% put = cnfdamput(Smax, dS, T, dT, X, R, SIG);
% Crank-Nicolson Finite Differences method for Pricing an American Put
% Smax : maximum stock price
% dS : increment of stock price
% T : maturity date
% dT : time step
% X : exercise price
% R : risk free interest rate
%
% reference : John C. Hull, Options, Futures, and Other Derivatives
% 3rd Ed., Chap 15
M = ceil(Smax/dS); ds = Smax / M;
Svec = 0:ds:Smax;
N = ceil(T/dT); dt = T / N;
Tvec = 0:dt:T;
J = 1:M-1;
a = (-.25*R*dt*J + .25*SIG^2*dt*J.^2);
b = (1 - 0.5*SIG^2*dt*J.^2 - .5*R*dt);
c = (.25*R*dt*J + .25*SIG^2*dt*J.^2);
A = diag(b) + diag(a(2:M-1), -1) + diag(c(1:M-2), 1);
aa = .25*R*dt*J - .25*SIG^2*dt*J.^2;
bb = 1 + .5*SIG^2*dt*J.^2 + .5*R*dt;
cc = -.25*R*dt*J - .25*SIG^2*dt*J.^2;
B = diag(bb) + diag(aa(2:M-1), -1) + diag(cc(1:M-2), 1);
g = zeros(1, M-1);
put = zeros(N+1, M+1);
put(N+1, :) = max(X - Svec, 0);
put(:, 1) = X;
put(:, M+1) = 0;
for i = N:-1:1
g = put(i+1, 2:M) * A'; 
g(1) = g(1) + a(1)*X - aa(1)*X;
put(i, 2:M) = [B \ g']';
put(i, :) = max(max(X - Svec,0), put(i,:));
end
mesh(Svec, Tvec,put);
end