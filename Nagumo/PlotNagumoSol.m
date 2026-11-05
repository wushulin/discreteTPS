% --------------------------------------------- %
% Script to plot target and exact control for
% the Nagumo optimal control problem. Figure 5.4
% --------------------------------------------- %
clear; clc; close all

gam = 5;

L = 20;
T = 5;

hx = 2^-2;
ht = 2^-3*hx^2;
t = 0:ht:T;
Nt = length(t) - 1;

%space discretization finite differences
[A, W, x] = SpaceFD(hx, L);
Nx = length(x);


%initial condition
y0 = 1.2*sqrt(3)*(x >= 9/20*L & x <= 11/20*L);


%target function
Ynat = SolveStateFD(A, ht, Nt, gam, zeros(Nx, Nt+1), y0);
ns = round(0.5*T/ht) + 1;
Yd = [Ynat(:, 1:ns), repmat(Ynat(:, ns), 1, Nt+1-ns)];


mesh(x,t,Yd'); 
view(2);
colorbar;
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);
xlabel('$x$', 'Interpret', 'latex', 'fontsize', 21);
ylabel('$t$', 'Interpret', 'latex', 'fontsize', 21);
ylim([0 T]);


%saveas(gca, 'fig/Nagumo_gam_5_L_20_T_5_hx_2^-2_ht_2^-3hx^2_target_2D', 'epsc');



%exact control
ys = Ynat(:, ns);
Uex = zeros(Nx, Nt+1);
Uex(:, ns+1:end) = repmat(A*ys + gam/3*ys.^3 - ys, 1, Nt+1-ns);


mesh(x,t,Uex'); 
view(2);
colorbar;
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);
xlabel('$x$', 'Interpret', 'latex', 'fontsize', 21);
ylabel('$t$', 'Interpret', 'latex', 'fontsize', 21);
ylim([0 T]);


%saveas(gca, 'fig/Nagumo_gam_5_L_20_T_5_hx_2^-2_ht_2^-3hx^2_u_exact_2D', 'epsc');