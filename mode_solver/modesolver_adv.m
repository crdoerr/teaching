%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finite difference mode solver
% Chris Doerr 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% simulation variables
lambdac = 1.310;     % vacuum center wavelength
Nx = 300;      % 700
Ny = 300;      % 700
xlim = 1.5;
ylim = 1.0;
%% waveguide parameters
wg.w = 0.5;         % waveguide width
wg.h = 0.22;        % waveguide height
ncore = 3.48;       % index of core
nclad = 1.45;       % index of cladding
neff_guess = 2.5;   % initial guess for effective index
%% grid calculations
x = linspace(-xlim,xlim,Nx);
y = linspace(-ylim,ylim,Ny);
Dx = x(2) - x(1);
Dy = y(2) - y(1);
%% fill in er
ermat = ones(Nx,Ny)*nclad^2;
%---Si waveguide
ermat = addwgpiece(-wg.w/2,-wg.h,wg.w/2,0,ncore,ermat,Dx,Dy);
%% calculate modes
[neff1,Ex1,Ey1,Ez1,Hx1,Hy1,Hz1] = modesolve_v2(lambdac,x,y,ermat,neff_guess,4,'TE',1);
shg