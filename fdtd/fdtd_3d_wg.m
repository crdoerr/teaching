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
Dx2 = Dx/2; Dy2 = Dy/2; Dz2 = Dz/2;
Nx2 = Nx*2+1; Ny2 = Ny*2+1; Nz2 = Nz*2+1;
er2 = noxide^2*ones(Nx2,Ny2,Nz2);
%---first straight section
for s = 0:Dx2/5:l
    for w = -wg.w/2:Dx2/5:wg.w/2
        for h = -wg.h/2:Dz2/5:wg.h/2
            xp = round(s/Dx2);
            yp = round((R0+l+w)/Dy2);
            zp = round((h+zlim)/Dz2);
            if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1 && zp > 0 && zp < Nz2+1
                er2(xp,yp,zp) = nsi^2;
            end
        end
    end
end
%---curved section
for theta = pi/2:-pi/4000:0
    for R = R0-wg.w/2:Dx2/5:R0+wg.w/2
        for h = -wg.h/2:Dz2/5:wg.h/2
            xp = round((R*cos(theta) + l)/Dx2);
            yp = round((R*sin(theta) + l)/Dy2);
            zp = round((h+zlim)/Dz2);
            if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1 && zp > 0 && zp < Nz2+1
                er2(xp,yp,zp) = nsi^2;
            end
        end
    end
end
%---second straight section
for s = 0:Dy2/5:l
    for w = -wg.w/2:Dy2/5:wg.w/2
        for h = -wg.h/2:Dz2/5:wg.h/2
            xp = round((R0+l+w)/Dx2);
            yp = round(s/Dy2);
            zp = round((h+zlim)/Dz2);
            if xp > 0 && xp < Nx2+1 && yp > 0 && yp < Ny2+1 && zp > 0 && zp < Nz2+1
                er2(xp,yp,zp) = nsi^2;
            end
        end
    end
end
for xn = 1:Nx;
    for yn = 1:Ny;
        for zn = 1:Nz;
            erx(xn,yn,zn) = ...
                (er2(2*xn,2*yn-1,2*zn-1) + er2(2*xn+1,2*yn-1,2*zn-1)...
                + er2(2*xn,2*yn,2*zn-1) + er2(2*xn+1,2*yn,2*zn-1)...
                + er2(2*xn,2*yn-1,2*zn) + er2(2*xn+1,2*yn-1,2*zn)...
                + er2(2*xn,2*yn,2*zn) + er2(2*xn+1,2*yn,2*zn))/8;
            ery(xn,yn,zn) = ...
                (er2(2*xn-1,2*yn,2*zn-1) + er2(2*xn-1,2*yn+1,2*zn-1)...
                + er2(2*xn,2*yn,2*zn-1) + er2(2*xn,2*yn+1,2*zn-1)...
                + er2(2*xn-1,2*yn,2*zn) + er2(2*xn-1,2*yn+1,2*zn)...
                + er2(2*xn,2*yn,2*zn) + er2(2*xn,2*yn+1,2*zn))/8;
            erz(xn,yn,zn) = ...
                (er2(2*xn-1,2*yn-1,2*zn) + er2(2*xn-1,2*yn-1,2*zn+1)...
                + er2(2*xn,2*yn-1,2*zn) + er2(2*xn,2*yn-1,2*zn+1)...
                + er2(2*xn-1,2*yn,2*zn) + er2(2*xn-1,2*yn,2*zn+1)...
                + er2(2*xn,2*yn,2*zn) + er2(2*xn,2*yn,2*zn+1))/8;
        end
    end
end
erplotxy = -abs(gradient(reshape(erx(:,:,round(Nz/2)),Nx,Ny)).'...
    + gradient(reshape(ery(:,:,round(Nz/2)),Nx,Ny).'));
erplotyz = -abs(gradient(reshape(ery(round(Nx/4),:,:),Ny,Nz)).'...
    + gradient(reshape(erz(round(Nx/4),:,:),Ny,Nz).'));
%% define profile of launch field
for yn = 1:Ny;
    for zn = 1:Nz;
        envyz(1,yn,zn) = 2100*Dx*exp(-((yn*Dy-R0-l)/wg.w*2)^2)...
            *exp(-((zn*Dz-zlim)/wg.h*2)^2);
    end
end
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
    Ey(2,:,:) = Ey(2,:,:) + envyz*cos(2*pi*fc*t);   %*exp(-((t-pulsew*3)/pulsew)^2);
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
        h = surf(x,y,Hz2D.' + erplotxy.');
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
        Hz2D = reshape(Hz(round(Nx/4),:,:),Ny,Nz);
        h = surf(y,z,Hz2D.' + erplotyz);
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