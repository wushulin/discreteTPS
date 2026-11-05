% --------------------------------------------- %
% Script for Newton outer loop and TPS inner loop
% to run each pair (nu, T, tht), with nu=1e-3,
% 1e-4, 1e-5, 1e-6, T=2^0, 2^3, 2^5, and
% tht = 0, 0.5, 1. One needs to change these 
% values manually to test each case. The results
% are stored in NagumoNewtonTPSdata.txt. Table 5.3
% --------------------------------------------- %
clear; clc; close all

nu = 1e-3;
gam = 5;

L = 20;
T = 2^0;

hx = 2^-2;
ht = 2^-3*hx^2;
t = 0:ht:T;
Nt = length(t) - 1;

tht = 0;
tol = 1e-10;


%space discretization finite differences
[A, W, x] = SpaceFD(hx, L);
Nx = length(x);


%initial condition
y0 = 1.2*sqrt(3)*(x >= 9/20*L & x <= 11/20*L);


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



%Newton+TPS
twl = 2^-1;
Ntw = T/twl
Nl = Nt/Ntw;

NitTPS = 1e1; NitNew = 5e1;

gl = zeros(Nx, Ntw-1); gr = gl;
Yk = zeros(Nx, Nt+1); Pk = Yk;
res = [];

for k = 1:NitNew
  k
  for nit = 1:NitTPS
    nit
    gltp = zeros(size(gl)); grtp = gltp;
    ytp = []; ptp = [];

    for m = 1:Ntw
      if m == 1
        [Ak1, Ck1, fy1, fp1] = SpaceNewtonFD(hx, L, gam, Yk(:, 1:Nl+1), Pk(:, 1:Nl+1));
        M1 = BuildThtMatrixNewton(Ak1, Ck1, ht, Nl, nu, tht, speye(Nx), sparse(Nx, Nx), sparse(Nx, Nx), speye(Nx));
        F1 = BuildThtRHSNewton(ht, Nl, tht, Yd(:, 1:Nl+1)', fy1, fp1, y0, gr(:, 1));
        X1 = M1\F1; X1 = reshape(X1, Nx, 2*Nl+2);
        y1 = X1(:, 1:Nl+1); p1 = X1(:, Nl+2:end);

        ytp = [ytp y1]; ptp = [ptp p1(:, 1:end-1)];
        gltp(:, 1) = y1(:, end);

      elseif m == Ntw
        [Akn, Ckn, fyn, fpn] = SpaceNewtonFD(hx, L, gam, Yk(:, end-Nl:end), Pk(:, end-Nl:end));
        Mn = BuildThtMatrixNewton(Akn, Ckn, ht, Nl, nu, tht, speye(Nx), sparse(Nx, Nx), sparse(Nx, Nx), speye(Nx));
        Fn = BuildThtRHSNewton(ht, Nl, tht, Yd(:, end-Nl:end)', fyn, fpn, gl(:, end), zeros(Nx, 1));      
        Xn = Mn\Fn; Xn = reshape(Xn, Nx, 2*Nl+2);
        yn = Xn(:, 1:Nl+1); pn = Xn(:, Nl+2:end);

        ytp = [ytp yn(:, 2:end)]; ptp = [ptp pn]; 
        grtp(:, end) = pn(:, 1);

      else
        [Akm, Ckm, fym, fpm] = SpaceNewtonFD(hx, L, gam, Yk(:, (m-1)*Nl+1:m*Nl+1), Pk(:, (m-1)*Nl+1:m*Nl+1));
        Mm = BuildThtMatrixNewton(Akm, Ckm, ht, Nl, nu, tht, speye(Nx), sparse(Nx, Nx), sparse(Nx, Nx), speye(Nx));
        Fm = BuildThtRHSNewton(ht, Nl, tht, Yd(:, (m-1)*Nl+1:m*Nl+1)', fym, fpm, gl(:, m-1), gr(:, m));
        Xm = Mm\Fm; Xm = reshape(Xm, Nx, 2*Nl+2);
        ym = Xm(:, 1:Nl+1); pm = Xm(:, Nl+2:end);
  
        ytp = [ytp ym(:, 2:end)]; ptp = [ptp pm(:, 1:end-1)];
        gltp(:, m) = ym(:, end); grtp(:, m-1) = pm(:, 1);
      end
    end

    gl = gltp; gr = grtp;

    if 0
    figure(1)
    mesh(t, x, ytp); 
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    xlabel('$t$', 'interpret', 'latex', 'fontsize', 21);
    ylabel('$x$', 'interpret', 'latex', 'fontsize', 21);

    figure(2)
    mesh(t, x, ptp/nu);
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    xlabel('$t$', 'interpret', 'latex', 'fontsize', 21);
    ylabel('$x$', 'interpret', 'latex', 'fontsize', 21);

    pause(0.5)
    end
  end

  res = [res max(abs([ytp-Yk; ptp-Pk]), [], 'all')]
  Yk = ytp; Pk = ptp;

  if res(end)<tol, break; end
end
YnewTPS = Yk; PnewTPS = Pk; UnewTPS = Pk/nu;

[sqrt(ip(YnewTPS-Yd, YnewTPS-Yd)) sqrt(ip(UnewTPS-Uex, UnewTPS-Uex))]

sqrt(ip(YnewTPS-Ynum, YnewTPS-Ynum)) 

sqrt(ip(UnewTPS-Unum, UnewTPS-Unum))