# Discrete analysis for time parallel Schwarz

This repository contains the data and the code for the paper *New time domain decomposition methods for parabolic optimal control problems III: Discrete analysis*.

## Main folder

Explanations for each `.m` file or `.zip` file are as follows: 

`Fig3_1_ConvV.m`----the code for plotting Figure 3.1, the condition number of the eigenvector matrix of the time discrete matrix

`Fig4_1_rho_and_its_bound_real_case_wrt_z.zip`----the codes for plotting Figure 4.1, that is the spectral radius of the iteration 
matrix R of the TPS method and its bound in the real eigenvalue case studied in Section 4.3.1. These codes are used to plot 
$\rho(R)$ and its bound for fixed penalty parameter $\nu$ and as $z=\Delta t\lambda$ varies. 

`Fig4_2_rho_and_its_bound_real_case_wrt_nu.zip`----similar to Figure 4.1, these codes are used to plot $\max_z\rho(R)$ and its bound 
$\max_z\sqrt{|\tau_1|^2+|\tau_2|^2}$ as the penalty parameter $\nu$ varies from 1e-7 to 1.

`Fig4_3 (bound_of_rho_complex_case).zip`----these codes are used for the contour plotting in Figure 4.3 in the case that the 
eigenvalues of the coefficient matrix of the state equation are complex. We use the bound $|\tau_1|+|\tau_2|$ analyzed in Section
4.3.2 for this plotting

`rho_and_its_bound_complex_case.zip`----these codes are not used to plot the figure appearing in the paper, but you can use them to plot $\rho(R)$ and its bound for a given penalty parameter $\nu$ 


## L-shape
This folder contains numerical tests for the optimal control problems constrained by an advection-diffusion problem:

$$\begin{aligned}
\partial_t y + \bm c \cdot \nabla y - \kappa \Delta y&=u& \text{in }& \Omega\times(0, T),\\
y(\bm x, t) &= 0& \text{on }& \partial\Omega\times(0, T),\\
y(\bm x, 0) &= y_0& \text{in }& \Omega.
\end{aligned}$$

`BuildThtMatrix2d.m`: function to build one big space-time matrix $M$ for the reduced optimality system $M\boldsymbol U=\bm F$, with $\boldsymbol U = [\boldsymbol Y, \boldsymbol P]$.

`BuildThtRHS2d.m`: function to build one big space-time right-hand side $\boldsymbol F$ for the reduced optimality system.

`LshapeTPSAdvDiffPar.m`: script to run the TPS for the advection-diffusion case for three $\theta$-methods and six $\nu$ values. It was run in parallel on the supercomputer LUNARC in Lund University and data was stored in the folder *data*. 

`LshapeTPSDiffPar.m`: script to run the TPS for the pure diffusion case for three $\theta$-methods and six $\nu$ values. It was run in parallel on the supercomputer LUNARC in Lund University and data was stored in the folder *data*. 

`PlotLshapeSol.m`: script to plot the target $\hat y$, the control $u$ and the solution $y$ on the L-shape domain.

`PlotLshapeTPSAdvDiff_Nbitvsnu.m`: script to plot the weak scalability test for the advection-diffusion case using the data generated from `LshapeTPSAdvDiffPar.m`.

`PlotLshapeTPSDiff_3dt.m`: script to plot the explicit Euler using three different time step $\Delta t$ in the pure diffusion case with two different values of $\nu$.

`PlotLshapeTPSDiff_3tht.m`: script to plot the three $\theta$-methods in the pure diffusion case for $\nu=$ 1e-6.

`PlotLshapeTPSDiff_Nbitvsnu.m`: script to plot the weak scalability test for the pure diffusion case using the data generated from `LshapeTPSDiffPar.m`.

`SpaceDisLshape.m`: function to build space matrix using a P1 finite element method.

`yhat2dmh.m`: function to build a target as a moving hot spot.

## Nagumo
This folder contains numerical tests for the optimal control problems constrained by a reaction-diffusion problem:

$$\begin{aligned}
\partial_t y - \partial_{xx} y + \frac{\gamma}3 y^3 - y &= u& \text{ in }& \Omega\times(0,T),\\
\partial_x y(0, t) = \partial_xy(L, t) &= 0& \text{ in }& (0, T),\\
y(x, 0) &= y_0(x)& \text{ in }& \Omega.
\end{aligned}$$

`BuildThtMatrixNewton.m`: function to build one big space-time matrix $M$ for one Newton step of the linearized reduced optimality system $M\boldsymbol U=\boldsymbol F$, with $\boldsymbol U = [\boldsymbol Y, \boldsymbol P]$.

`BuildThtRHSNewton.m`: function to build one big space-time right-hand side $\bm F$ for one Newton step of the linearized reduced optimality system.

`NagumoNewtonTPS.m`: script for Newton outer loop and TPS inner loop to test for three $\theta$-method, three number of subdomains $N$, and four values of $\nu$. The data are saved in `NagumoNewtonTPSdata.txt`

`NagumoTPSNewton.m`: script for TPS outer loop and Newton inner loop to test for three $\theta$-method, three number of subdomains $N$, and four values of $\nu$. One needs to run this script, which can be run in parallel, to get the data and save in a seperate folder `data`. We cannot upload them due to the size.

`PlotNagumoSol.m`: script to plot the target $\hat y$ and the exact control $u$, i.e., $\nu=0$.

`PlotNagumoTPSNewton.m`: script to plot the weak scalability test using the data generated from `NagumoTPSNewton.m`. 

`SolveStateFD.m`: function to solve the forward reaction-diffusion problem using finite difference in space and explicit Euler in time. 

`SpaceFD.m`: function to generate the Laplace matrix with Neumann boundary condition. 

`SpaceNewtonFD.m`: function to generate space matrix for one Newton step. 