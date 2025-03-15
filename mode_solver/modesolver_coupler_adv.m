%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finite difference mode solver
% Chris Doerr 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% simulation variables
lambdac = 1.310;     % vacuum center wavelength
Nx = 200;      % 700
Ny = 200;      % 700
xlim = 2.0;
ylim = 1.0;
%% waveguide parameters
wg.w = 0.4;         % waveguide width
wg.h = 0.22;        % waveguide height
wg.sp = 1.2;        % center-to-center spacing
ncore = 3.48;       % index of core
nclad = 1.45;       % index of cladding
neff_guess = 2.9;   % initial guess for effective index
%% grid calculations
x = linspace(-xlim,xlim,Nx);
y = linspace(-ylim,ylim,Ny);
Dx = x(2) - x(1);
Dy = y(2) - y(1);
%% fill in er
ermat = ones(Nx,Ny)*nclad^2;
%---Si waveguides
ermat = addwgpiece(-wg.w/2-wg.sp/2,-wg.h,wg.w/2-wg.sp/2,0,ncore,ermat,Dx,Dy);
ermat = addwgpiece(-wg.w/2+wg.sp/2,-wg.h,wg.w/2+wg.sp/2,0,ncore,ermat,Dx,Dy);
%% calculate modes
[neff1,Ex1,Ey1,Ez1,Hx1,Hy1,Hz1] = modesolve_v2(lambdac,x,y,ermat,neff_guess,4,'TE',1);
shg