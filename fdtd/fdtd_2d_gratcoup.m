%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-D FDTD simulator
% (no variation in z direction)
% Chris Doerr
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
nsi = 3.50;           % index of silicon
%% adjustable variables
lambdac = 1.310;       % vacuum center wavelength
time = 200e-15;       % time length of simulation (s)
xlim = 7;             % limits of x window (window width = 2*xlim)
ylim = 5;             % limits of y window (window width = 2*ylim)
Dx = 0.02;            % x grid size
Dy = Dx;              % y grid size
pulsew = 10e-15;      % input pulse width (s)
Nplot = 20;            % plots every Nplot time steps (to speed execution)
%% waveguide parameters
wg.h = 0.34;            % waveguide height
grat.pitch = 0.51;      % grating pitch
grat.depth = 0.120;      % grating depth
grat.num = 20;          % number of grating notches
grat.start = 2;         % length of waveguide before starting grating
l = xlim*2;             % length of straight waveguide
%% calculated variables
Dt = Dx/(c0/noxide)/1.5;
x = -xlim:Dx:xlim;
y = -ylim:Dy:ylim;
Nx = length(x);
Ny = length(y);
numsteps = round(time/Dt);
fc = c0/lambdac;
%% set up structure
Dx2 = Dx/2; Dy2 = Dy/2;
Nx2 = Nx*2+1; Ny2 = Ny*2+1;
er2 = noxide^2*ones(Nx2,Ny2);
%---straight waveguide
for s = 0:Dx2/5:l
    for h = -wg.h/2:Dx2/5:wg.h/2
        xp = round(s/Dx2);
        yp = round((2.5+wg.h/2+h)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
    end
end
%---grating
for m = 1:grat.num
    for s = -grat.pitch/4:Dx2/5:grat.pitch/4
        for h = -grat.depth:Dy2/5:0
            xp = round((s+grat.start+(m-1)*grat.pitch)/Dx2);
            yp = round((2.5+wg.h+h)/Dy2);
            if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
                er2(xp,yp) = noxide^2;
            end
        end
    end
end
%% calculate erx, ery
for xn = 1:Nx
    for yn = 1:Ny
        erx(xn,yn) = (er2(2*xn,2*yn-1) + er2(2*xn+1,2*yn-1)...
            + er2(2*xn,2*yn) + er2(2*xn+1,2*yn))/4;
        ery(xn,yn) = (er2(2*xn-1,2*yn) + er2(2*xn-1,2*yn+1)...
            + er2(2*xn,2*yn) + er2(2*xn,2*yn+1))/4;
        erz(xn,yn) = ...
            (er2(2*xn-1,2*yn-1) ...
            + er2(2*xn,2*yn-1) ...
            + er2(2*xn-1,2*yn) ...
            + er2(2*xn,2*yn))/4;
    end
end
erplot = -abs(gradient(erx).' + gradient(ery.'));
%% define profile of launch field
ya = (1:Ny)*Dy;
envy = 50*Dx*exp(-((ya-2.5-wg.h/2)/wg.h*2/2).^2)/2;
%% initialize fields
Ex = zeros(Nx,Ny);
Ey = zeros(Nx,Ny);
Ez = zeros(Nx,Ny);
Hx = zeros(Nx,Ny);
Hy = zeros(Nx,Ny);
Hz = zeros(Nx,Ny);
%% calculate constant combinations to speed execution
Dtu0Dy = Dt/u0/Dy;
Dtu0Dx = Dt/u0/Dx;
Dte0x = Dt./(e0*erx);
Dte0y = Dt./(e0*ery);
Dte0z = Dt./(e0*erz);
%% main loop
fprintf('Simulation starting...  ')
for cnt = 1:numsteps
    t = cnt*Dt;
    %% inject field
    Ez(2,1:Ny) = Ez(2,1:Ny) + envy.*cos(2*pi*fc*t);  %*exp(-((t-pulsew*3)/pulsew)^2);
    %% calculate H field
    Hxnew = Hx  - Dtu0Dy*(circshift(Ez,[0 -1]) - Ez);
    Hynew = Hy + Dtu0Dx*(circshift(Ez,[-1 0]) - Ez);
    Hynew(Nx,:) = -Ez(Nx,:)./sqrt(u0/e0./erz(Nx,:));
    Hxnew(:,Ny) = Ez(:,Ny)./sqrt(u0/e0./erz(:,Ny));
    Hx = Hxnew;
    Hy = Hynew;
    %% calculate E field
    Eznew = Ez + Dte0z.*((Hy - circshift(Hy,[1 0]))/Dx - (Hx - circshift(Hx,[0 1]))/Dy);
    Eznew(1,:) = Hy(1,:).*sqrt(u0/e0./ery(1,:));
    Eznew(:,1) = -Hx(:,1).*sqrt(u0/e0./erx(:,1));
    Ez = Eznew;
    %% plot field
    if round(cnt/Nplot)==cnt/Nplot
        h = surf(x,y,2*Ez.'+ erplot);
        axis([-xlim xlim -ylim ylim -inf inf])
        axis equal
        clim([-1 1])
        view(0,90)
        set(h,'linestyle','none')
        set(gca,'fontsize',14)
        title('\fontsize{14}E_z')
        xlabel('\fontsize{14}x (\mum)')
        ylabel('\fontsize{14}y (\mum)')
        colormap jet
        drawnow
        shg
    end
end
fprintf('Simulation finished.\n')
shg