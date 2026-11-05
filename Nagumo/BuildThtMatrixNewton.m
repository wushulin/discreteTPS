function M = BuildThtMatrixNewton(Ak, Ck, dt, N, nu, tht, Ay0, Ap0, AyT, ApT)
% Build the matrix M to solve the problem MU = F
% with U = [Y; P] of one Newton step for the 
% Nagumo optimality system.
% Space discretization: finite difference with
% the matrix Ak linearized operator and Ck 
% coupling in the adjoint part.
% Time discretization: linear theta-method:
% y_{n+1} - y_n + dt*(tht*A_{n+1}*y_{n+1} 
% + (1-tht)*A_n*y_n) - dt/nu*(tht*p_{n+1} 
% + (1-tht)*p_n) = dt*(tht*fy_{n+1} 
% + (1-tht)*fy_n),
% p_n - p_{n+1} + dt*(tht*A_n*p_n 
% + (1-tht)*A_{n+1}*p_{n+1}) + dt*(tht*C_n*y_n
% + (1-tht)*C_{n+1}*y_{n+1}) = dt*(tht*g_n 
% + (1-tht)*g_{n+1}),
% with g = yhat + fp, where A_n and C_n are the 
% blocks of Ak and Ck at the time point t_n 
% from SpaceNewtonFD.m 
% theta=0 for explicit Euler, theta=0.5 for 
% trapezoidal rule, and theta=1 for implicit 
% Euler.
% There are four components in the matrix
% M = [M11, M12; M21, M22]. The first block of 
% M11 and M12 and the last block of M21 and M22 
% are specially adapted to different transmission
% conditions. 
%
% Input
% ----------
% Ak : matrix
%     block diagonal linearized operator at 
%     the N+1 time points
% Ck : matrix
%     block diagonal coupling of y in the 
%     adjoint equation
% dt : float
%     time step
% N : int
%     number of points in time (N+1)
% nu : float
%     penalization for control
% tht : float
%     theta value
% Ay0 : matrix
%     a (J+1)*(J+1) matrix for first block of M11
% Ap0 : matrix
%     a (J+1)*(J+1) matrix for first block of M12
% AyT : matrix
%     a (J+1)*(J+1) matrix for last block of M21
% ApT : matrix
%     a (J+1)*(J+1) matrix for last block of M22

nk = size(Ak, 1);
nx = nk/(N+1);

e = ones(N+1, 1);
I = speye(nk);
S = kron(spdiags(e, -1, N+1, N+1), speye(nx));

%construct four matrices
M11 = I + tht*dt*Ak - S*(I - (1-tht)*dt*Ak);
M22 = I + tht*dt*Ak - S'*(I - (1-tht)*dt*Ak);
M12 = -dt/nu*(tht*I + (1-tht)*S);
M21 = dt*(tht*Ck + (1-tht)*S'*Ck);

%adapt the initial and final conditions
M11(1:nx, 1:nx) = Ay0;
M12(1:nx, 1:nx) = Ap0;
M21(end-nx+1:end, end-nx+1:end) = AyT;
M22(end-nx+1:end, end-nx+1:end) = ApT;

%concatenation of four matrices
M = [M11, M12; M21, M22];
end