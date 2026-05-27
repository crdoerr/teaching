%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finite difference mode solver
% Chris Doerr 2022; Updated Rishay Gupta 2026
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
ncore = 3.50;       % index of core
nclad = 1.45;       % index of cladding
neff_guess = 2.5;   % initial guess for effective index
%% grid calculations
x = unique([-xlim:0.04:-0.7, -0.7:0.01:0.7, 0.7:0.04:xlim]); % coarse spacing on the edges of the guide
y = unique([-ylim:0.04:-0.4, -0.4:0.01:0.4, 0.4:0.04:ylim]); % dense spacing around the center of the guide

Nx = length(x);
Ny = length(y);
%% fill in er
ermat = ones(Nx,Ny)*nclad^2;
%---Si waveguide
ermat = addwgpiece(-wg.w/2,-wg.h,wg.w/2,0,ncore,ermat,x,y);
%% calculate modes
[neff1,Ex1,Ey1,Ez1,Hx1,Hy1,Hz1] = modesolve_v2(lambdac,x,y,ermat,neff_guess,4,'TE',1);
