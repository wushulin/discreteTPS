function Y = SolveStateFD(A, dt, N, gam, U, y0)
% Solve the Nagumo equation 
% y_t - y_xx + (gam/3*y^3 - y) = u
% forward in time by explicit Euler with the 
% finite difference matrix A from SpaceFD.m
% y_{n+1} = y_n - dt*(A*y_n + f(y_n)) + dt*u_n,
% f(y) = gam/3*y.^3 - y.
% It is stable for dt <= 2/(4/h^2 + max f'(y)).
%
% Input
% ----------
% A : matrix
%     finite difference matrix
% dt : float
%     time step
% N : int
%     number of points in time (N+1)
% gam : float
%     parameter in the nonlinear reaction term
% U : matrix
%     control at the nodes (time along the columns)
% y0 : vector
%     initial condition
%
% Output
% -------
% Y : matrix
%     state at the nodes (time along the columns)

f = @(y) gam/3*y.^3 - y;

Y = zeros(length(y0), N+1);
Y(:, 1) = y0;

for n = 1:N
  y = Y(:, n);
  Y(:, n+1) = y - dt*(A*y + f(y)) + dt*U(:, n);
end
end