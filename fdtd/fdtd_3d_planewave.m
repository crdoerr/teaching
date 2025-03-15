%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-D FDTD simulator
% Chris Doerr 2010
% Bell Labs, Alcatel-Lucent
% all length units in um
%
%       ^ y
%       |
%       |
%      /----> x
%     /
%    L z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
clf
%% physical constants
e0 = 8.854188e-12*1e-6;
u0 = 4*pi*1e-7*1e-6;
c0 = 1/sqrt(e0*u0);
%% material parameters
noxide = 1.45;        % index of oxide
nsi = 3.48;           % index of silicon
%% adjustable variables
lambdac = 1.55;       % vacuum center wavelength
time = 200e-15;       % time length of simulation (s)
xlim = 4;             % limits of x window (window width = 2*xlim)
ylim = 4;             % limits of y window (window width = 2*ylim)
zlim = 0.5;           % limits of z window (window width = 2*zlim)
Dx = 0.05;            % x grid size
Dy = Dx;              % y grid size
Dz = Dx;              % z grid size
pulsew = 10e-15;      % input pulse width (s)
Nplot = 5;            % plots every Nplot time steps (to speed execution)
%% waveguide parameters
wg.w = 1;             % waveguide width
wg.h = 0.22;          % waveguide height
R0 = 4;               % bend radius
l = 2;                % length of straight waveguide
%% calculated variables
Dt = Dx/(c0/noxide)/1.9;
x = -xlim:Dx:xlim;
y = -ylim:Dy:ylim;
z = -zlim:Dz:zlim;
Nx = length(x);
Ny = length(y);
Nz = length(z);
numsteps = round(time/Dt);
fc = c0/lambdac;
%% set up structure
erx = ones(Nx,Ny,Nz)*noxide^2;
ery = ones(Nx,Ny,Nz)*noxide^2;
erz = ones(Nx,Ny,Nz)*noxide^2;
%% initialize fields
Ex = zeros(Nx,Ny,Nz);
Ey = zeros(Nx,Ny,Nz);
Ez = zeros(Nx,Ny,Nz);
Hx = zeros(Nx,Ny,Nz);
Hy = zeros(Nx,Ny,Nz);
Hz = zeros(Nx,Ny,Nz);
%% calculate constant combinations to speed execution
Dtu0Dx = Dt/u0/Dx;
Dtu0Dy = Dt/u0/Dy;
Dtu0Dz = Dt/u0/Dz;
Dte0x = Dt./(e0*erx);
Dte0y = Dt./(e0*ery);
Dte0z = Dt./(e0*erz);
%% main loop
fprintf('Simulation starting...  ')
colormap jet
for cnt = 1:numsteps
    t = cnt*Dt;
    %% inject field
    Hz(round(Nx/4),:,:) = Hz(round(Nx/4),:,:)...
        + 0.5*cos(2*pi*fc*t);
    Ey(round(Nx/4),:,:) = Ey(round(Nx/4),:,:)...
        + 0.5*sqrt(u0/(e0*noxide^2))*cos(2*pi*fc*t);
    %% calculate H field
    Hxnew = Hx + Dtu0Dz*(circshift(Ey,[0 0 -1]) - Ey) - Dtu0Dy*(circshift(Ez,[0 -1 0]) - Ez);
    Hynew = Hy + Dtu0Dx*(circshift(Ez,[-1 0 0]) - Ez) - Dtu0Dz*(circshift(Ex,[0 0 -1]) - Ex);
    Hznew = Hz + Dtu0Dy*(circshift(Ex,[0 -1 0]) - Ex) - Dtu0Dx*(circshift(Ey,[-1 0 0]) - Ey);
    Hznew(Nx,:,:) = Ey(Nx,:,:)./sqrt(u0/e0./ery(Nx,:,:));
    Hynew(Nx,:,:) = -Ez(Nx,:,:)./sqrt(u0/e0./erz(Nx,:,:));
    Hxnew(:,Ny,:) = Ez(:,Ny,:)./sqrt(u0/e0./erz(:,Ny,:));
    Hznew(:,Ny,:) = -Ex(:,Ny,:)./sqrt(u0/e0./erx(:,Ny,:));
    Hxnew(:,:,Nz) = -Ey(:,:,Nz)./sqrt(u0/e0./ery(:,:,Nz));
    Hynew(:,:,Nz) = Ex(:,:,Nz)./sqrt(u0/e0./erx(:,:,Nz));
    Hx = Hxnew;
    Hy = Hynew;
    Hz = Hznew;
    %% calculate E field
    Exnew = Ex + Dte0x.*((Hz - circshift(Hz,[0 1 0]))/Dy - (Hy - circshift(Hy,[0 0 1]))/Dz);
    Eynew = Ey + Dte0y.*((Hx - circshift(Hx,[0 0 1]))/Dz - (Hz - circshift(Hz,[1 0 0]))/Dx);
    Eznew = Ez + Dte0z.*((Hy - circshift(Hy,[1 0 0]))/Dx - (Hx - circshift(Hx,[0 1 0]))/Dy);
    Eznew(1,:,:) = Hy(1,:,:).*sqrt(u0/e0./ery(1,:,:));
    Eynew(1,:,:) = -Hz(1,:,:).*sqrt(u0/e0./erz(1,:,:));
    Exnew(:,1,:) = Hz(:,1,:).*sqrt(u0/e0./erz(:,1,:));
    Eznew(:,1,:) = -Hx(:,1,:).*sqrt(u0/e0./erx(:,1,:));
    Exnew(:,:,1) = -Hy(:,:,1).*sqrt(u0/e0./ery(:,:,1));
    Eynew(:,:,1) = Hx(:,:,1).*sqrt(u0/e0./erx(:,:,1));
    Ex = Exnew;
    Ey = Eynew;
    Ez = Eznew;
    %% plot field
    if round(cnt/Nplot)==cnt/Nplot
        subplot(121)
        Hz2D = reshape(Hz(:,:,round(Nz/2)),Nx,Ny);
        h = surf(x,y,Hz2D.');
        axis([-xlim xlim -ylim ylim -inf inf])
        axis equal
        caxis([-1 1])
        view(0,90)
        set(h,'linestyle','none')
        set(gca,'fontsize',14)
        title('\fontsize{14}H_z')
        xlabel('\fontsize{14}x (\mum)')
        ylabel('\fontsize{14}y (\mum)')
        subplot(122)
        Hz2D = reshape(Hz(round(Nx/2),:,:),Ny,Nz);
        h = surf(y,z,Hz2D.');
        axis([-xlim xlim -ylim ylim -inf inf])
        axis equal
        caxis([-1 1])
        view(0,90)
        set(h,'linestyle','none')
        set(gca,'fontsize',14)
        title('\fontsize{14}H_z')
        xlabel('\fontsize{14}y (\mum)')
        ylabel('\fontsize{14}z (\mum)')
        drawnow
        shg
    end
end
fprintf('Simulation finished.\n')
shg