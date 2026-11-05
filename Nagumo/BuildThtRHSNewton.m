function F = BuildThtRHSNewton(dt, N, tht, Yhat, fy, fp, F0, FT)
% Build the vector F to solve the problem MU = F
% of one Newton step for the Nagumo optimality
% system. The matrix M is generated from
% BuildThtMatrixNewton.m 
% There are two components in the vector 
% F = [F1; F2]. F1 corresponds to the right-hand 
% side of the state Y with the source fy; 
% F2 corresponds to the right-hand side of the 
% adjoint P with the target and the source fp.
%
% Input
% ----------
% dt : float
%     time step
% N : int
%     number of points in time (N+1)
% tht : float
%     theta value
% Yhat : matrix
%     (N+1)*(J+1) target function, time 
%     along the rows
% fy : matrix
%     (J+1)*(N+1) source of the state equation 
%     from SpaceNewtonFD.m
% fp : matrix
%     (J+1)*(N+1) source of the adjoint equation 
%     from SpaceNewtonFD.m
% F0 : vector
%     initial condition
% FT : vector
%     final condition

F1 = dt*(tht*fy(:, 2:end) + (1-tht)*fy(:, 1:end-1));
F1 = [F0(:); F1(:)];

G = Yhat' + fp;
F2 = dt*(tht*G(:, 1:end-1) + (1-tht)*G(:, 2:end));
F2 = [F2(:); FT(:)];

F = [F1; F2];
end