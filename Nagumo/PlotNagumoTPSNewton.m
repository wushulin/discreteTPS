% --------------------------------------------- %
% Script to generate the table (or plot) for the 
% weak scalability test for three time 
% integrators with different values of \nu in the
% Nagumo case. N varies in {2,16,64}.
% To change the plot, one needs to change in the
% title nu=1e-3, 1e-6 and in saveas nu=1e-3, 1e-6
% accordingly. Figure 5.5
% --------------------------------------------- %
clear; clc; close all

TT = 2.^[0 3 5];
Tht = (0:2)/2;
markers = {'o','s','h'};


figure('Position', [100, 100, 800, 600]);
for j = 1:length(TT)
  for tht = Tht
    filename = sprintf('data/NagumoTPSNewton_gam_5_L_20_hx_2^-2_ht_2^-3hx^2_twl_2^-1_NitTPS_20_tol_1e-10_nu_1e-6_Ntw_%d_tht_%g.mat', TT(j)/2^-1, tht);

    S = load(filename);
    err = (S.erry+S.errp)/(S.erry(1)+S.errp(1));

    semilogy(1:length(err), err, '-', 'Marker', markers{2*tht+1}, 'linewidth', 3, 'markersize', 20);
    hold on
  end
end
hold off
set(gca, 'fontsize', 20);
set(gca, 'linewidth', 2.5);
xlabel('Iter', 'interpret', 'latex', 'fontsize', 25);
xlim([1 10]);
%ylim([1e-2 1e1]); %in the case nu=1e-3
ylabel('Error', 'interpret', 'latex', 'fontsize', 25);
title('$\nu = 10^{-6}$', 'interpret', 'latex', 'fontsize', 21);
legend({'$\theta = 0, N = 2$', '$\theta = 0.5, N = 2 $', '$\theta = 1, N = 2 $', '$\theta = 0, N = 16$', '$\theta = 0.5, N = 16$', '$\theta = 1, N = 16$', '$\theta = 0, N = 64$', '$\theta = 0.5, N = 16$', '$\theta = 1, N = 64$'}, 'interpreter', 'latex', 'location', 'best', 'fontsize', 21);
grid on
ax = gca;
ax.YMinorGrid = 'off';


%saveas(gca, 'fig/NagumoTPSNewton_gam_5_L_20_hx_2^-2_ht_2^-3hx^2_twl_2^-1_NitTPS_10_tol_1e-10_nu_1e-6_3Ntw3tht_9curves', 'epsc');