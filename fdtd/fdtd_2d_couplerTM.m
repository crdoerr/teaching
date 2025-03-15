%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-D FDTD simulator
% (no variation in z direction)
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
nsi = 2.06;           % TM effective index of silicon
%% adjustable variables
lambdac = 1.55;       % vacuum center wavelength
time = 2000e-15;       % time length of simulation (s)
xlim = 19;             % limits of x window (window width = 2*xlim)
ylim = 3;             % limits of y window (window width = 2*ylim)
Dx = 0.05;            % x grid size
Dy = Dx;              % y grid size
pulsew = 10e-15;      % input pulse width (s)
Nplot = 10;            % plots every Nplot time steps (to speed execution)
%% waveguide parameters
wg.w = 0.5;             % waveguide width
wg.gap = 0.5;
R0 = 10;               % bend radius
l = 28;                % length of straight waveguide
a = 15*pi/180;
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
offset = wg.w/2 + wg.gap/2 + 2*R0*(1-cos(a));
%---first bend
for theta = pi/2:-pi/4000:pi/2-a
    for R = R0-wg.w/2:Dx2/5:R0+wg.w/2
        xp = round((R*cos(theta))/Dx2);
        yp = round((R*sin(theta)+(offset-R0)+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
        yp = round((-R*sin(theta)-(offset-R0)+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
    end
end
%---second bend
for theta = -pi/2-a:pi/4000:-pi/2
    for R = R0-wg.w/2:Dx2/5:R0+wg.w/2
        xp = round((R*cos(theta)+R0*2*sin(a))/Dx2);
        yp = round((R*sin(theta)+(offset-2*R0*(1-cos(a))+R0)+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
        yp = round((-R*sin(theta)-(offset-2*R0*(1-cos(a))+R0)+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
    end
end
%---straight sections
for s = 0:Dx2/5:l
    for w = -wg.w/2:Dx2/5:wg.w/2
        xp = round((s+2*R0*sin(a))/Dx2);
        yp = round((wg.gap/2+wg.w/2+w+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
        yp = round((-wg.gap/2-wg.w/2-w+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
    end
end
%---third bend
for theta = -pi/2:pi/4000:-pi/2+a
    for R = R0-wg.w/2:Dx2/5:R0+wg.w/2
        xp = round((R*cos(theta)+R0*2*sin(a)+l)/Dx2);
        yp = round((R*sin(theta)+(offset-2*R0*(1-cos(a))+R0)+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
        yp = round((-R*sin(theta)-(offset-2*R0*(1-cos(a))+R0)+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
    end
end
%---fourth bend
for theta = pi/2+a:-pi/4000:pi/2
    for R = R0-wg.w/2:Dx2/5:R0+wg.w/2
        xp = round((R*cos(theta)+l+4*R0*sin(a))/Dx2);
        yp = round((R*sin(theta)+(offset-R0)+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
        yp = round((-R*sin(theta)-(offset-R0)+ylim)/Dy2);
        if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1
            er2(xp,yp) = nsi^2;
        end
    end
end
%% calculate erx, ery, erz
for xn = 1:Nx;
    for yn = 1:Ny;
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
envy = 50*Dx*exp(-((ya-offset-ylim)/wg.w*2/2).^2)/2;
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
colormap jet
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
        h = surf(x,y,Ez.'+ erplot);
        axis([-xlim xlim -ylim ylim -inf inf])
        axis equal
        caxis([-1 1])
        view(0,90)
        set(h,'linestyle','none')
        set(gca,'fontsize',14)
        title('\fontsize{14}E_z')
        xlabel('\fontsize{14}x (\mum)')
        ylabel('\fontsize{14}y (\mum)')
        drawnow
        shg
    end
end
fprintf('Simulation finished.\n')
shg