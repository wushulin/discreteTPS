function [Ak, Ck, fy, fp, W, x] = SpaceNewtonFD(h, L, gam, yk, pk)
% Finite difference space matrices of one 
% Newton step for the Nagumo optimality system 
% on (0, L) with homogeneous Neumann BC,
% linearized at (yk, pk):
%   y_t + Ak*y - p/nu = fy,
%  -p_t + Ak*p + Ck*y = yhat + fp,
% with Ak = A + diag(gam*yk.^2 - 1), 
% Ck = diag(1 + 2*gam*yk.*pk),
% fy = 2*gam/3*yk.^3, fp = 2*gam*yk.^2.*pk 
% and A, W, x from SpaceFD.m
% yk, pk are given at one time point, (J+1)*1, 
% or at all time points, (J+1)*(N+1); then Ak 
% and Ck are block diagonal over the time points
% (time-major) and fy, fp have same shape as yk.
%
% Input
% ----------
% h : float
%     space meshsize (L/h must be an integer)
% L : float
%     length of the interval
% gam : float
%     parameter in the nonlinear reaction term
% yk : vector or matrix
%     state of the previous Newton iterate 
%     at the nodes
% pk : vector or matrix
%     adjoint of the previous Newton iterate 
%     at the nodes
%
% Output
% -------
% Ak : matrix
%     linearized operator A'(yk)
% Ck : matrix
%     coupling of y in the adjoint equation, 
%     1 + A''(yk)*pk
% fy : vector or matrix
%     source of the state equation A'(yk)*yk-A(yk)
% fp : vector or matrix
%     source of the adjoint equation A''(yk)*pk*yk
% W : matrix
%     trapezoidal weights
% x : vector
%     nodes

[A, W, x] = SpaceFD(h, L);
Nn = size(yk, 2);
nk = numel(yk);

Ak = kron(speye(Nn), A) + spdiags(gam*yk(:).^2 - 1, 0, nk, nk);
Ck = spdiags(1 + 2*gam*yk(:).*pk(:), 0, nk, nk);
fy = 2*gam/3*yk.^3;
fp = 2*gam*yk.^2.*pk;
end