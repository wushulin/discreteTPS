function [Ms, K, C, msh] = SpaceDisLshape(h, b)
% Generate a uniform P1 finite element
% discretization of the L-shape domain 
% Omega = (-1,1)^2 \ [0,1]x[-1,0] with meshsize h,
% and assemble the finite element matrices of the
% advection-diffusion operator -kap*Delta y +
% b.grad y with homogeneous Dirichlet BCs.
% Each square of size h is split into two 
% triangles by its diagonal.
% Ms, K and C are restricted to the interior 
% nodes, so the spatial operator is A = kap*K + C,
% and the semi-discrete state equation reads
% Ms*y' + A*y = Ms*u.
%
% Input
% ----------
% h : float
%     space meshsize, 1/h must be an integer
% c : function handle or vector
%     advection field c(x, y) = [c1, c2] with x, y column vectors and
%     a two-column output, or a constant 1x2 vector (default [0 0])
%
% Output
% -------
% Ms : matrix
%     Nin*Nin mass matrix (phi_j, phi_i)
% K : matrix
%     Nin*Nin stiffness matrix 
%     (grad phi_j, grad phi_i)
% C : matrix
%     Nin*Nin advection matrix 
%     (b.grad phi_j, phi_i)
% msh : struct
%     p   : Np*2 node coordinates
%     tri : Ne*3 triangles (counterclockwise)
%     in  : Nin*1 indices of interior nodes 
%           (the unknowns)
%     bd  : indices of boundary nodes
%     Mf  : Nin*Np mass matrix (phi_j, phi_i) 
%           for interior i and all j, to build 
%           load vectors from nodal values f on 
%.          all nodes: Mf*f
%     h   : meshsize

if nargin < 2
  b = [0 0];
end

if isnumeric(b)
  c0 = b(:)';
  c = @(x, y) repmat(c0, numel(x), 1);
end

%integer grid on (-1,1)^2, m is the index of the reentrant corner
n = round(2/h);
m = n/2;
if abs(n*h-2) > 1e-12 || mod(n, 2) ~= 0
  error('1/h must be an integer.');
end
[ix, iy] = ndgrid(0:n, 0:n);

%remove the nodes in the quadrant (0,1]x[-1,0)
keep = ~(ix > m & iy < m);
Np = nnz(keep);
ID = zeros(n+1, n+1);
ID(keep) = 1:Np;
p = [-1+ix(keep)*h, -1+iy(keep)*h];

%boundary nodes: outer square and the two reentrant edges
ix = ix(keep); iy = iy(keep);
isbd = ix == 0 | ix == n | iy == 0 | iy == n ...
       | (ix == m & iy <= m) | (iy == m & ix >= m);
bd = find(isbd);
in = find(~isbd);

%triangles: split each kept square along the diagonal n1-n3
[sx, sy] = ndgrid(0:n-1, 0:n-1);
sq = ~(sx >= m & sy < m);
sx = sx(sq); sy = sy(sq);
n1 = ID(sub2ind([n+1 n+1], sx+1, sy+1));
n2 = ID(sub2ind([n+1 n+1], sx+2, sy+1));
n3 = ID(sub2ind([n+1 n+1], sx+2, sy+2));
n4 = ID(sub2ind([n+1 n+1], sx+1, sy+2));
tri = [n1 n2 n3; n1 n3 n4];
Ne = size(tri, 1);

%area and gradients of the barycentric functions
x1 = p(tri(:, 1), 1); x2 = p(tri(:, 2), 1); x3 = p(tri(:, 3), 1);
y1 = p(tri(:, 1), 2); y2 = p(tri(:, 2), 2); y3 = p(tri(:, 3), 2);
ar = ((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1))/2;
gx = [y2-y3, y3-y1, y1-y2]./(2*ar);
gy = [x3-x2, x1-x3, x2-x1]./(2*ar);

%advection field at the edge midpoints, midpoint q is opposite vertex q
xm = [(x2+x3)/2, (x3+x1)/2, (x1+x2)/2];
ym = [(y2+y3)/2, (y3+y1)/2, (y1+y2)/2];
bq = b(xm(:), ym(:));
bx = reshape(bq(:, 1), Ne, 3);
by = reshape(bq(:, 2), Ne, 3);

%local matrices, the advection uses the edge
%midpoint rule (exact for linear b), where 
%phi_i = 1/2 at the midpoints q ~= i and 0 
%at q = i
I = zeros(Ne, 9); J = I; VM = I; VK = I; VC = I;
k = 0;
for i = 1:3
  for j = 1:3
    k = k + 1;
    I(:, k) = tri(:, i); J(:, k) = tri(:, j);
    VM(:, k) = ar/12*(1 + (i == j));
    VK(:, k) = ar.*(gx(:, i).*gx(:, j) + gy(:, i).*gy(:, j));
    q = setdiff(1:3, i);
    VC(:, k) = ar/6.*sum(bx(:, q).*gx(:, j) + by(:, q).*gy(:, j), 2);
  end
end

%global assembly
Ms = sparse(I(:), J(:), VM(:), Np, Np);
K = sparse(I(:), J(:), VK(:), Np, Np);
C = sparse(I(:), J(:), VC(:), Np, Np);

%mesh structure
msh.p = p;
msh.tri = tri;
msh.in = in;
msh.bd = bd;
msh.Mf = Ms(in, :);
msh.h = h;

%restrict to the interior nodes
Ms = Ms(in, in);
K = K(in, in);
C = C(in, in);
end