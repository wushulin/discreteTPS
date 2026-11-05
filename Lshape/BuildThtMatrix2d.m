function M = BuildThtMatrix2d(Ms, A, dt, N, nu, al, tht, Ay0, Ap0, AyT, ApT)
% Build the matrix M to solve the forward-backward 
% linear system MU = F with U = [Y; P]. 
% Space discretization: P1 finite element on the 
% L-shape domain with the mass matrix Ms and the 
% spatial matrix A = kap*K + C, all generated from 
% SpaceDisLshape.m
% Time discretization: linear theta-method, where % theta=0 for explicit Euler, theta=0.5 for 
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
% Ms : matrix
%     a Nin*Nin mass matrix
% A : matrix
%     a Nin*Nin spatial operator
% dt : float
%     time step
% N : int
%     number of points in time (N+1)
% nu : float
%     penalization for control
% al : float
%     penalization for target
% tht : float
%     theta value: 0, 0.5, 1
% Ay0 : matrix
%     a Nin*Nin matrix for the first block of M11
% Ap0 : matrix
%     a Nin*Nin matrix for the first block of M12
% AyT : matrix
%     a Nin*Nin matrix for the last block of M21
% ApT : matrix
%     a Nin*Nin matrix for the last block of M22

%number of unknowns in space
Nin = size(Ms, 1);

%construct four matrices
e = ones(N+1, 1);
M11 = kron(speye(N+1), Ms+dt*tht*A) ...
      -kron(spdiags(e, -1, N+1, N+1), Ms-(1-tht)*dt*A);
M22 = M11';
M21 = kron(speye(N+1), tht*dt*Ms) ...
      +kron(spdiags(e, 1, N+1, N+1), (1-tht)*dt*Ms);
M12 = -M21'/nu;
M21 = al*M21;

%adapt the initial and final conditions
M11(1:Nin, 1:Nin) = Ay0;
M12(1:Nin, 1:Nin) = Ap0;
M21(end-Nin+1:end, end-Nin+1:end) = AyT;
M22(end-Nin+1:end, end-Nin+1:end) = ApT;

%concatenation of four matrices
M = [M11, M12; M21, M22];
end