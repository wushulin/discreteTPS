% --------------------------------------------- %
% Script to plot the target, solution and control
% for moving hot spot in an L-shape domain.
% Figures 5.1 and 5.2
% --------------------------------------------- %
clear; clc; close all

nu = 1e-1
al = 1;
%kap = 2^-1; b = [0 0]; %case I
kap = 2^-3; b = [1 1]; %case II

hx = 2^-4;
T = 2^1;
ht = 0.5*hx;
t = 0:ht:T;
Nt = length(t) - 1;

tht = 1;


%spatial discretization
[Ms, K, C, msh] = SpaceDisLshape(hx, b);
A = kap*K + C;
Nin = length(msh.in);
Np = size(msh.p, 1);

%target function
Yhat = yhat2dmh(msh.p(:, 1)', msh.p(:, 2)', t', kap, b, T);

%solve entire system
M = BuildThtMatrix2d(Ms, A, ht, Nt, nu, al, tht, speye(Nin), sparse(Nin, Nin), sparse(Nin, Nin), speye(Nin));
dM = decomposition(M, 'lu');
F = BuildThtRHS2d(msh.Mf, ht, Nt, tht, al*Yhat, Yhat(1, msh.in)', sparse(Nin,1));
X = dM\F; X = reshape(X, Nin, 2*Nt+2);
y = X(:, 1:Nt+1); p = X(:, Nt+2:end);
u = p/nu;


%select specific time to plot solution
%[~, n] = min(abs(t - 0*T));
%[~, n] = min(abs(t - 0.5*T));
[~, n] = min(abs(t - 1*T));
yf = zeros(Np, 1); yf(msh.in) = y(:, n);
uf = zeros(Np, 1); uf(msh.in) = u(:, n);


%%y
trisurf(msh.tri, msh.p(:, 1), msh.p(:, 2), yf);
shading interp; 
view(2); 
axis equal tight; 
colorbar;
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);
title([sprintf('$y(%g T)$', t(n)/T)], 'Interpret', 'latex', 'fontsize', 21);
xlabel('$x_1$', 'Interpret', 'latex', 'fontsize', 21);
ylabel('$x_2$', 'Interpret', 'latex', 'fontsize', 21);


%saveas(gca, 'fig/LshapeMovehotspot_heat_y_05T', 'epsc');

%saveas(gca, 'fig/LshapeMovehotspot_advdiff_y_1T', 'epsc');


%%u
trisurf(msh.tri, msh.p(:, 1), msh.p(:, 2), uf);
shading interp; 
view(2); 
axis equal tight; 
colorbar;
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);
title([sprintf('$u(%g T)$', t(n)/T)], 'Interpret', 'latex', 'fontsize', 21);
xlabel('$x_1$', 'Interpret', 'latex', 'fontsize', 21);
ylabel('$x_2$', 'Interpret', 'latex', 'fontsize', 21);
%clim([0 1]); %use when t=T


%saveas(gca, 'fig/LshapeMovehotspot_heat_u_05T', 'epsc');

%saveas(gca, 'fig/LshapeMovehotspot_advdiff_u_1T', 'epsc');


%%target
trisurf(msh.tri, msh.p(:, 1), msh.p(:, 2), Yhat(n,:));
shading interp; 
view(2); 
axis equal tight; 
colorbar;
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);
title([sprintf('$\\hat y(%g T)$', t(n)/T)], 'Interpret', 'latex', 'fontsize', 21);
xlabel('$x_1$', 'Interpret', 'latex', 'fontsize', 21);
ylabel('$x_2$', 'Interpret', 'latex', 'fontsize', 21);


%saveas(gca, 'fig/LshapeMovehotspot_heat_target_05T', 'epsc');

%saveas(gca, 'fig/LshapeMovehotspot_advdiff_target_1T', 'epsc');