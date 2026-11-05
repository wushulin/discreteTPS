% --------------------------------------------- %
% Script to plot the figure for three time 
% integrators to compare the convergence behavior
% in the pure diffusion case. Figure 5.3 left
% --------------------------------------------- %
clear; clc; close all

nu = 1e-6;
al = 1;
kap = 2^-1;

hx = 2^-3;
T = 2^0;
ht = 2^-3*hx^2;
t = 0:ht:T;
Nt = length(t) - 1;
	
Twl = 2^-1;
NTwl = T/Twl;

Tht = (0:2)/2;
Nit = 2e1;

erry = zeros(Nit, length(Tht));
errp = erry;

%spatial discretization
[Ms, K, C, msh] = SpaceDisLshape(hx);
A = kap*K + C;
Nin = length(msh.in);
Np = size(msh.p, 1);

rng(0);
for tht = Tht
  tht

  %local matrix
  Mm = BuildThtMatrix2d(Ms, A, ht, Nt/NTwl, nu, al, tht, speye(Nin), sparse(Nin, Nin), sparse(Nin, Nin), speye(Nin));
  dMm = decomposition(Mm, 'lu');
  
  %initiate transmission condition
  gl = rand(Nin, NTwl-1); gr = gl;

  for nit = 1:Nit
    nit

    %allocate temporal values
    gltp = zeros(size(gl)); grtp = gltp;
    ytp = []; ptp = []; 

    for m = 1:NTwl
      if m == 1

        %solve first time window
        F1 = BuildThtRHS2d(msh.Mf, ht, Nt/NTwl, tht, al*sparse(Nt/NTwl+1, Np), sparse(Nin, 1), gr(:, 1));
        X1 = dMm\F1; X1 = reshape(X1, Nin, 2*Nt/NTwl+2);
        y1 = X1(:, 1:Nt/NTwl+1); p1 = X1(:, Nt/NTwl+2:end);

        %save solution and transmission data
        ytp = [ytp y1]; ptp = [ptp, p1(:, 1:end-1)];
        gltp(:, 1) = y1(:, end);

      elseif m == NTwl

        %solve last time window
        Fn = BuildThtRHS2d(msh.Mf, ht, Nt/NTwl, tht, al*sparse(Nt/NTwl+1, Np), gl(:, end), sparse(Nin, 1));
        Xn = dMm\Fn; Xn = reshape(Xn, Nin, 2*Nt/NTwl+2);
        yn = Xn(:, 1:Nt/NTwl+1); pn = Xn(:, Nt/NTwl+2:end);

        %save solution and transmission data
        ytp = [ytp yn(:,2:end)]; ptp = [ptp, pn];
        grtp(:, end) = pn(:, 1);

      else 

        %solve inner time window
        Fm = BuildThtRHS2d(msh.Mf, ht, Nt/NTwl, tht, al*sparse(Nt/NTwl+1, Np), gl(:, m-1), gr(:, m));
        Xm = dMm\Fm; Xm = reshape(Xm, Nin, 2*Nt/NTwl+2);
        ym = Xm(:, 1:Nt/NTwl+1); pm = Xm(:, Nt/NTwl+2:end);

        %save solution and transmission data
        ytp = [ytp ym(:,2:end)]; ptp = [ptp, pm(:, 1:end-1)];
        gltp(:, m) = ym(:, end); grtp(:, m-1) = pm(:, 1);
      end
    end

    %update transmission conditions
    gl = gltp; gr = grtp;

    erry(nit, tht*2+1) = norm(ytp, 'inf');
    errp(nit, tht*2+1) = norm(ptp, 'inf');
  end
end

semilogy(1:2:Nit, (erry(1:2:end, 1)+errp(1:2:end, 1))/(erry(1, 1)+errp(1, 1)), '--o', 1:2:Nit, (erry(1:2:end, 2)+errp(1:2:end, 2))/(erry(1, 2)+errp(1, 2)), '--s', 1:2:Nit, (erry(1:2:end, 3)+errp(1:2:end, 3))/(erry(1, 3)+errp(1, 3)), '--h', 'MarkerSize', 12, 'linewidth', 2);
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);
xlabel('Iteration', 'interpret', 'latex', 'fontsize', 21);
ylabel('Error', 'interpret', 'latex', 'fontsize', 21);
legend({'$\theta = 0$', '$\theta = 0.5$', '$\theta = 1$'}, 'interpreter', 'latex', 'location', 'best', 'fontsize', 21);
xlim([1 Nit]);
grid on
ax = gca;
ax.YMinorGrid = 'off';


saveas(gca, 'fig/LshapeTPSDiff_nu_1e-6_kap_2^-1_hx_2^-3_ht_2^-3hx^2_T_2^0_Twl_2^-1_Nit_20_3tht', 'epsc');