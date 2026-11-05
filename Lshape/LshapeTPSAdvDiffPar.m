% --------------------------------------------- %
% Script to run on the supercomputer LUNARC 
% for a Slurm job array and generate data for
% PlotLshapeTPSAdvDiff_Nbitvsnu.m
% Each array task solves one (nu, tht) pair,
% chosen from SLURM_ARRAY_TASK_ID = 1, ..., 18,
% and saves its own file.
% Only run this script for T=2^5 with Nit=180, 
% or remove the test queue option in the bash
% file, because LUNARC cannot finish with Nit=200
% using the test queue option. No problem for
% T=2^0 and 2^3. Table 5.2
% --------------------------------------------- %
clear; clc; close all

al = 1;
kap = 2^-3;
b = [1, 1];

hx = 2^-3;
T = 2^5;
ht = 2^-3*hx^2;
t = 0:ht:T;
Nt = length(t) - 1;

Twl = 2^-1;
NTwl = T/Twl

Nit = 180;%2e2;
Tht = (0:2)/2;
Nu = 10.^-(1:1:6);


%choose one (nu, tht) pair from the Slurm array index 1, ..., 18
idx = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(idx), idx = 1; end
[i, j] = ind2sub([length(Tht) length(Nu)], idx);
tht = Tht(i)
nu = Nu(j)

erry = zeros(Nit, 1);
errp = erry;

%spatial discretization
[Ms, K, C, msh] = SpaceDisLshape(hx, b);
A = kap*K + C;
Nin = length(msh.in);
Np = size(msh.p, 1);

%local matrix
Mm = BuildThtMatrix2d(Ms, A, ht, Nt/NTwl, nu, al, tht, speye(Nin), sparse(Nin, Nin), sparse(Nin, Nin), speye(Nin));
dMm = decomposition(Mm, 'lu');

%initiate transmission condition, same random start for all pairs
rng(0);
gl = rand(Nin, NTwl-1); gr = gl;

for nit = 1:Nit
  nit
  gltp = zeros(size(gl)); grtp = gltp;
  ytp = []; ptp = [];

  for m = 1:NTwl
    if m == 1
      F1 = BuildThtRHS2d(msh.Mf, ht, Nt/NTwl, tht, al*sparse(Nt/NTwl+1, Np), sparse(Nin, 1), gr(:, 1));
      X1 = dMm\F1; X1 = reshape(X1, Nin, 2*Nt/NTwl+2);
      y1 = X1(:, 1:Nt/NTwl+1); p1 = X1(:, Nt/NTwl+2:end);
      ytp = [ytp y1]; ptp = [ptp, p1(:, 1:end-1)];
      gltp(:, 1) = y1(:, end);

    elseif m == NTwl
      Fn = BuildThtRHS2d(msh.Mf, ht, Nt/NTwl, tht, al*sparse(Nt/NTwl+1, Np), gl(:, end), sparse(Nin, 1));
      Xn = dMm\Fn; Xn = reshape(Xn, Nin, 2*Nt/NTwl+2);
      yn = Xn(:, 1:Nt/NTwl+1); pn = Xn(:, Nt/NTwl+2:end);
      ytp = [ytp yn(:,2:end)]; ptp = [ptp, pn];
      grtp(:, end) = pn(:, 1);

    else
      Fm = BuildThtRHS2d(msh.Mf, ht, Nt/NTwl, tht, al*sparse(Nt/NTwl+1, Np), gl(:, m-1), gr(:, m));
      Xm = dMm\Fm; Xm = reshape(Xm, Nin, 2*Nt/NTwl+2);
      ym = Xm(:, 1:Nt/NTwl+1); pm = Xm(:, Nt/NTwl+2:end);
      ytp = [ytp ym(:,2:end)]; ptp = [ptp, pm(:, 1:end-1)];
      gltp(:, m) = ym(:, end); grtp(:, m-1) = pm(:, 1);
    end
  end

  gl = gltp; gr = grtp;
  erry(nit) = norm(ytp, 'inf');
  errp(nit) = norm(ptp, 'inf');
end

%save one file for this (nu, tht) pair
filename = sprintf('LshapeTPSAdvDiff_kap_2^-3_b_1_hx_2^-3_ht_2^-3hx^2_T_2^5_Twl_2^-1_Nit_%d_nu_1e-%d_tht_%g.mat', Nit, j, tht)
save(filename, 'al', 'kap', 'hx', 'T', 'ht', 'Twl', 'NTwl', 'Nit', 'nu', 'tht', 'erry', 'errp');
fprintf("Save completed. \n");