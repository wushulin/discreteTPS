function F = BuildThtRHS2d(Mf, dt, N, tht, Yhat, F0, FT)
% Build the vector F to solve the problem MU = F 
% with U = [Y; P]. The matrix M is defined in
% BuildThtMatrix2d.m
% There are two components in the vector 
% F = [F1; F2]. F1 corresponds to the right-hand 
% side of the state Y, and F2 corresponds to the
% right-hand side of the adjoint P. 
% The target is weighted by the mass matrix Mf.
%
% Input
% ----------
% Mf : matrix
%     a Nin*Np mass matrix for the target, 
%     e.g. msh.Mf with Yhat on all nodes, or Ms 
%     with Yhat on the interior nodes (Np = Nin)
% dt : float
%     time step
% N : int
%     number of points in time (N+1)
% tht : float
%     theta value
% Yhat : matrix
%     target function \hat y(x, t) of size 
%    (N+1)*Np, time along the rows
% F0 : vector
%     initial condition
% FT : vector
%     final condition

%number of unknowns in space
Nin = size(Mf, 1);

%construct F1 to solve Y
F1 = zeros((N+1)*Nin, 1);

%add the initial condition
F1(1:Nin) = F0;

%construct F2 to solve P
F2 = dt*Mf*((1-tht)*Yhat(2:end, :) + tht*Yhat(1:end-1, :))';

%add the final condition
F2 = [F2(:); FT(:)];

%concatenate F1 and F2
F = [F1; F2];
end