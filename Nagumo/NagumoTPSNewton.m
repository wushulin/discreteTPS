% --------------------------------------------- %
% Script for TPS outer loop and Newton inner loop
% to run on the supercomputer LUNARC for a Slurm
% job array and generate data for 
% PlotNagumoTPSNewton.m
% Each array task solves one (nu, T, tht) pair,
% chosen from SLURM_ARRAY_TASK_ID = 1, ..., 36,
% and saves its own file. Figure 5.5
% --------------------------------------------- %
clear; clc; close all

gam = 5;

L = 20;
hx = 2^-2;
ht = 2^-3*hx^2;

tol = 1e-10;

%space discretization finite differences
[A, W, x] = SpaceFD(hx, L);
Nx = length(x);

%initial condition
y0 = 1.2*sqrt(3)*(x >= 9/20*L & x <= 11/20*L);

Nu = 10.^-(3:6);
TT = 2.^[0 3 5];
Tht = (0:2)/2;

%choose one (nu, T, tht) pair from the Slurm array index 1, ..., 36
idx = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(idx), idx = 1; end
[i, j, l] = ind2sub([length(Nu) length(TT) length(Tht)], idx);
nu = Nu(i);
T = TT(j);
tht = Tht(l);

t = 0:ht:T;
Nt = length(t) - 1;

%target function
Ynat = SolveStateFD(A, ht, Nt, gam, zeros(Nx, Nt+1), y0);
ns = round(0.5*T/ht) + 1;
Yd = [Ynat(:, 1:ns), repmat(Ynat(:, ns), 1, Nt+1-ns)];

%exact control
ys = Ynat(:, ns);
Uex = zeros(Nx, Nt+1);
Uex(:, ns+1:end) = repmat(A*ys + gam/3*ys.^3 - ys, 1, Nt+1-ns);

%weighted L^2 inner product
w = ones(1, Nt+1); w(1) = 1-tht; w(end) = tht;
ip = @(a, b) ht*sum(w.*sum(a.*(W*b), 1));


%Newton for global solution
NitNew = 10; Yk = zeros(Nx, Nt+1); Pk = Yk;
for k = 1:NitNew
  [Ak, Ck, fy, fp] = SpaceNewtonFD(hx, L, gam, Yk, Pk);
  M = BuildThtMatrixNewton(Ak, Ck, ht, Nt, nu, tht, speye(Nx), sparse(Nx, Nx), sparse(Nx, Nx), speye(Nx));
  F = BuildThtRHSNewton(ht, Nt, tht, Yd', fy, fp, y0, zeros(Nx, 1));
  X = M\F; X = reshape(X, Nx, 2*Nt+2);
  Ytp = X(:, 1:Nt+1); Ptp = X(:, Nt+2:end);

  res = max(abs([Ytp-Yk; Ptp-Pk]), [], 'all')
  Yk = Ytp; Pk = Ptp;
  
  if res<tol, break; end
end
Ynum = Yk; Pnum = Pk; Unum = Pk/nu;

[sqrt(ip(Ynum-Yd, Ynum-Yd)) sqrt(ip(Unum-Uex, Unum-Uex))]


%TPS+Newton
NitTPS = 2e1;
NitNew = 2e1;
twl = 2^-1;
Ntw = T/twl
Nl = Nt/Ntw;

erry = ones(1, NitTPS); errp = erry;
gl = zeros(Nx, Ntw-1); gr = gl;
Yk = zeros(Nx, Nt+1); Pk = Yk;

for nit = 1:NitTPS
  gltp = zeros(size(gl)); grtp = gltp;
  ytp = []; ptp = [];

  for m = 1:Ntw
    if m == 1
      yk = Yk(:, 1:Nl+1); pk = Pk(:, 1:Nl+1);

      for k = 1:NitNew
        [Ak1, Ck1, fy1, fp1] = SpaceNewtonFD(hx, L, gam, yk, pk);
        M1 = BuildThtMatrixNewton(Ak1, Ck1, ht, Nl, nu, tht, speye(Nx), sparse(Nx, Nx), sparse(Nx, Nx), speye(Nx));
        F1 = BuildThtRHSNewton(ht, Nl, tht, Yd(:, 1:Nl+1)', fy1, fp1, y0, gr(:, 1));
        X1 = M1\F1; X1 = reshape(X1, Nx, 2*Nl+2);
        y1 = X1(:, 1:Nl+1); p1 = X1(:, Nl+2:end);

        res1 = max(abs([yk-y1; pk-p1]), [], 'all')
        yk = y1; pk = p1;

        if res1<tol, break; end
      end

      ytp = [ytp y1]; ptp = [ptp p1(:, 1:end-1)];
      gltp(:, 1) = y1(:, end);

    elseif m == Ntw
      yk = Yk(:, end-Nl:end); pk = Pk(:, end-Nl:end);

      for k = 1:NitNew
        [Akn, Ckn, fyn, fpn] = SpaceNewtonFD(hx, L, gam, yk, pk);
        Mn = BuildThtMatrixNewton(Akn, Ckn, ht, Nl, nu, tht, speye(Nx), sparse(Nx, Nx), sparse(Nx, Nx), speye(Nx));
        Fn = BuildThtRHSNewton(ht, Nl, tht, Yd(:, end-Nl:end)', fyn, fpn, gl(:, end), zeros(Nx, 1));      
        Xn = Mn\Fn; Xn = reshape(Xn, Nx, 2*Nl+2);
        yn = Xn(:, 1:Nl+1); pn = Xn(:, Nl+2:end);

        resn = max(abs([yk-yn; pk-pn]), [], 'all')
        yk = yn; pk = pn;

        if resn<tol, break; end
      end

      ytp = [ytp yn(:, 2:end)]; ptp = [ptp pn];  
      grtp(:, end) = pn(:, 1);

    else
      yk = Yk(:, (m-1)*Nl+1:m*Nl+1); pk = Pk(:, (m-1)*Nl+1:m*Nl+1);

      for k = 1:NitNew
        [Akm, Ckm, fym, fpm] = SpaceNewtonFD(hx, L, gam, yk, pk);
        Mm = BuildThtMatrixNewton(Akm, Ckm, ht, Nl, nu, tht, speye(Nx), sparse(Nx, Nx), sparse(Nx, Nx), speye(Nx));
        Fm = BuildThtRHSNewton(ht, Nl, tht, Yd(:, (m-1)*Nl+1:m*Nl+1)', fym, fpm, gl(:, m-1), gr(:, m));
        Xm = Mm\Fm; Xm = reshape(Xm, Nx, 2*Nl+2);
        ym = Xm(:, 1:Nl+1); pm = Xm(:, Nl+2:end);

        resm = max(abs([yk-ym; pk-pm]), [], 'all')
        yk = ym; pk = pm;

        if resm<tol, break; end
      end

      ytp = [ytp ym(:, 2:end)]; ptp = [ptp pm(:, 1:end-1)];
      gltp(:, m) = ym(:, end); grtp(:, m-1) = pm(:, 1);
    end
  end

  gl = gltp; gr = grtp;
  Yk = ytp; Pk = ptp;

  erry(1, nit) = sqrt(ip(ytp-Ynum, ytp-Ynum));
  errp(1, nit) = sqrt(ip(ptp-Pnum, ptp-Pnum));
end


%save one file for this (nu, tht) pair
filename = sprintf('NagumoTPSNewton_gam_5_L_20_hx_2^-2_ht_2^-3hx^2_twl_2^-1_NitTPS_20_tol_1e-10_nu_1e-%d_Ntw_%d_tht_%g.mat', i+2, Ntw, tht)
save(filename, 'gam', 'L', 'hx', 'tol', 'ht', 'T', 'nu', 'tht', 'twl', 'Ntw', 'NitTPS', 'NitNew', 'Ynum', 'Pnum', 'Unum', 'ytp', 'ptp', 'erry', 'errp');
fprintf("Save completed. \n");