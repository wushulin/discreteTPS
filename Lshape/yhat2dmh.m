function g = yhat2dmh(x, y, t, kap, c, T)
% Target function as a moving hot (mh) spot on 
% the L-shape domain. It starts as
% 10*exp(-50*|x - x0|^2) at x0 = (-0.5, 0) and is
% transported with the constant velocity 
% c = [c1 c2] and spreads by the diffusion kap. 
% The advection-diffusion solution in the whole 
% plane 10*s0/s(r)*exp(-|x - x0 - c*r|^2/s(r)),
% with s(r) = s0 + 4*kap*r and s0 = 1/50. 
% The time is rescaled r = t/T so that for any 
% final time T, the target at t = T is the same 
% as the one for T = 1. 
% The target is multiplied by a cutoff that is 
% 1 at distance >= bt from the boundary and 
% 0 on the boundary, so the target satisfies 
% the homogeneous Dirichlet boundary condition.
% Called with row vectors x, y of the node
% coordinates and a column vector t. 
% It returns a (Nt+1)*Np matrix, time along the
% rows, as used in BuildThtRHS2d.m.
%
% Input
% ----------
% x, y : vector
%     node coordinates
% t : vector
%     time points
% kap : float
%     diffusion coefficient
% c : vector
%     constant advection field [c1 c2]
% T : float
%     final time (default 1)

if nargin < 6
  T = 1;
end

x0 = [-0.5 0];
s0 = 1/50;
bt = 0.2;

%moving and spreading hot spot in the rescaled time 
r = t/T;
s = s0 + 4*kap*r;
g = 10*s0./s.*exp(-((x - x0(1) - c(1)*r).^2 + (y - x0(2) - c(2)*r).^2)./s);

%distance to the boundary of the L-shape (vertices counterclockwise)
V = [-1 -1; 0 -1; 0 0; 1 0; 1 1; -1 1; -1 -1];
d = inf(size(x));
for k = 1:6
  ex = V(k+1, 1) - V(k, 1); ey = V(k+1, 2) - V(k, 2);
  q = min(max(((x - V(k, 1))*ex + (y - V(k, 2))*ey)/(ex^2 + ey^2), 0), 1);
  d = min(d, sqrt((x - V(k, 1) - q*ex).^2 + (y - V(k, 2) - q*ey).^2));
end

%smooth cutoff: 0 on the boundary and 1 at distance >= bt
B = min(d/bt, 1);
g = g.*(B.^2.*(3 - 2*B));
end