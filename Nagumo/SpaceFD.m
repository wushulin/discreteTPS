function [A, W, x] = SpaceFD(h, L)
% Second-order finite difference approximation 
% of -\partial_xx on (0, L) with homogeneous
% Neumann BC by ghost points, at the J+1 nodes 
% x_j = j*h including the two end points: 
% row j is (-y_{j-1} + 2*y_j - y_{j+1})/h^2,
% and the two end rows are (2*y_0 - 2*y_1)/h^2 
% and (2*y_J - 2*y_{J-1})/h^2.
% W = h*diag(1/2, 1, ..., 1, 1/2) are the 
% trapezoidal weights, so that
% W*A is symmetric and the discrete L2 inner
% product is (a, b) = a'*W*b.
%
% Input
% ----------
% h : float
%     space meshsize (L/h must be an integer)
% L : float
%     length of the interval
%
% Output
% -------
% A : matrix
%     (J+1)*(J+1) finite difference matrix
% W : matrix
%     (J+1)*(J+1) diagonal matrix of 
%     trapezoidal weights
% x : vector
%     (J+1)*1 nodes

J = round(L/h);
if abs(J*h-L) > 1e-12*L
  error('L/h must be an integer.');
end

x = (0:J)'*h;
e = ones(J+1, 1);

A = spdiags([-e 2*e -e], -1:1, J+1, J+1)/h^2;
A(1, 2) = -2/h^2; A(end, end-1) = -2/h^2;

W = spdiags(h*e, 0, J+1, J+1);
W(1, 1) = h/2; W(end, end) = h/2;
end