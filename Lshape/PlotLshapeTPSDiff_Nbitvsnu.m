% --------------------------------------------- %
% Script to generate the table (or plot) for the 
% weak scalability test for three time 
% integrators with different values of \nu in the
% pure diffusion case. N varies in {2,16,64}.
% To change the plot, one needs to change in the 
% filename T_2^0,3,5, in title N=2,16,64, and 
% in saveas T_2^0,3,5 accordingly. Table 5.1
% --------------------------------------------- %
clear; clc; close all

Nu = 10.^-(1:1:6);
Tht = (0:2)/2;
markers = {'o','s','h'};

for tht = Tht
  Nind = [];
  for j = 1:6
    filename = sprintf('data/LshapeTPSDiff_kap_2^-1_hx_2^-3_ht_2^-3hx^2_T_2^0_Twl_2^-1_Nit_100_nu_1e-%d_tht_%g.mat',j,tht);

    S = load(filename);
    err = (S.erry+S.errp)/(S.erry(1)+S.errp(1));
    Nind = [Nind find(err<1e-4, 1)];
  end

  semilogx(Nu(1:length(Nind)),Nind,'-', 'Marker', markers{2*tht+1}, 'linewidth', 2, 'markersize', 15);
  hold on
end
hold off
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);
xlabel('Penalization parameter ($\nu$)', 'interpret', 'latex', 'fontsize', 21);
ylabel('Number of iteration', 'interpret', 'latex', 'fontsize', 21);
title('$N = 2$', 'interpret', 'latex', 'fontsize', 21);
legend({'$\theta = 0$', '$\theta = 0.5$', '$\theta = 1$'}, 'interpreter', 'latex', 'location', 'best', 'fontsize', 21);
xlim([1e-6 1e-1]);
grid on
ax = gca;
ax.XMinorGrid = 'off';


%saveas(gca, 'fig/LshapeTPSDiff_kap_2^-1_hx_2^-3_ht_2^-3hx^2_T_2^0_Twl_2^-1_Nit_100_Nbitvsnu', 'epsc');