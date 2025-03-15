function [neff_s,Exmat_s,Eymat_s,Ezmat_s,Hxmat_s,Hymat_s,Hzmat_s] = modesolve_v2(lambdac,x_,y_,ermat,neff_guess,nummodes,pol,modenum)
% function that calculates waveguide eigenmode
%% physical constants
e0 = 8.854188e-12*1e-6;
u0 = 4*pi*1e-7*1e-6;
c0 = 1/sqrt(e0*u0);
Z0 = sqrt(u0/e0);
%% calculated variables
Nx = length(x_);
Ny = length(y_);
Dx = x_(2) - x_(1);
Dy = y_(2) - y_(1);
k0 = 2*pi/lambdac;
%% calculate erx, ery, erz
%---improve in the future so that ermat has 2x the resolution
ermatx = (ermat + circshift(ermat,[-1,0]))/2;
ermaty = (ermat + circshift(ermat,[0,-1]))/2;
ermatz = ermat;
for y = 1:Ny
    erdiagx((y-1)*Nx+1:(y)*Nx) = ermatx(1:Nx,y);
    erdiagy((y-1)*Nx+1:(y)*Nx) = ermaty(1:Nx,y);
    erdiagz((y-1)*Nx+1:(y)*Nx) = ermatz(1:Nx,y);
end
erx = spdiags(erdiagx.',0,Nx*Ny,Nx*Ny);
ery = spdiags(erdiagy.',0,Nx*Ny,Nx*Ny);
erz = spdiags(erdiagz.',0,Nx*Ny,Nx*Ny);
%% set up matrices
temp = ones(Nx*Ny,1);
Ux = spdiags([-temp temp],[0 1],Nx*Ny,Nx*Ny)/Dx;
Uy = spdiags([-temp temp],[0 Nx],Nx*Ny,Nx*Ny)/Dy;
Uy(end,end) = 1/Dy;
P1 = [Ux/erz*Uy.' k0^2*speye(Nx*Ny)-Ux/erz*Ux.'; -k0^2*speye(Nx*Ny)+Uy/erz*Uy.' -Uy/erz*Ux.'];
P2 = [-Ux.'*Uy -k0^2*ery+Ux.'*Ux; k0^2*erx-Uy.'*Uy Uy.'*Ux];
P = P1*P2/k0^4;
%% find eigenvalues
options.disp = 1;
options.tol = 1e-9;
[V,D] = eigs(P,nummodes,neff_guess^2,options);
%% sort the eigenmodes from highest eigenvalue to lowest
[neff_,sorti] = sort(sqrt(diag(D)),'descend');
Vsort = 0*V;
for i = 1:length(neff_)
    Vsort(:,i) = V(:,sorti(i));
end
V = Vsort;
%% find desired mode
cnt = 0;
i = 0;
for m = 1:length(neff_)
    Ex = (V(1:Nx*Ny,m));
    Ey = (V(Nx*Ny+1:2*Nx*Ny,m));
    for y = 1:Ny
        Exmat(:,y) = real(Ex((y-1)*Nx+1:y*Nx));
        Eymat(:,y) = real(Ey((y-1)*Nx+1:y*Nx));
    end
    if sum(sum(abs(Exmat))) + sum(sum(abs(Eymat)))...
            < 4*(sum(sum(abs(Exmat(round(Nx/4):round(3*Nx/4),round(Ny/4):round(3*Ny/4)))))...
            + sum(sum(abs(Eymat(round(Nx/4):round(3*Nx/4),round(Ny/4):round(3*Ny/4))))))     % 2*
        if sum(sum(abs(Ex))) > sum(sum(abs(Ey))) && strcmp(pol,'TE')
            cnt = cnt + 1;
            if cnt==modenum
                i = m;
            end
        end
        if sum(sum(abs(Ex))) < sum(sum(abs(Ey))) && strcmp(pol,'TM')
            cnt = cnt + 1;
            if cnt==modenum
                i = m;
            end
        end
    end
end
%% calculate all field components (assuming x-y axes)
i = modenum;
neff = neff_(i);
Ex = (V(1:Nx*Ny,i));
Ey = (V(Nx*Ny+1:2*Nx*Ny,i));
Hz = 1i/k0/Z0*(-Uy*Ex + Ux*Ey);
Hy = erx*Ex/neff/Z0 - 1i*lambdac/2/pi/neff*Uy.'*Hz;
Hx = -ery*Ey/neff/Z0 - 1i*lambdac/2/pi/neff*Ux.'*Hz;
Ez = -1i/k0*inv(erz)*Z0*(-Ux.'*Hy + Uy.'*Hx);
%% convert to matrices
maxfield = max(max(abs(Ex)),max(abs(Ey)))*0.5;
Exmat_s = double(zeros(Nx,Ny));
Eymat_s = double(zeros(Nx,Ny));
Ezmat_s = double(zeros(Nx,Ny));
Hxmat_s = double(zeros(Nx,Ny));
Hymat_s = double(zeros(Nx,Ny));
Hzmat_s = double(zeros(Nx,Ny));
for y = 1:Ny
    Exmat_s(:,y) = real(Ex((y-1)*Nx+1:y*Nx))/maxfield;
    Eymat_s(:,y) = real(Ey((y-1)*Nx+1:y*Nx))/maxfield;
    Ezmat_s(:,y) = imag(Ez((y-1)*Nx+1:y*Nx))/maxfield;
    Hxmat_s(:,y) = real(Hx((y-1)*Nx+1:y*Nx))/maxfield;
    Hymat_s(:,y) = real(Hy((y-1)*Nx+1:y*Nx))/maxfield;
    Hzmat_s(:,y) = imag(Hz((y-1)*Nx+1:y*Nx))/maxfield;
end
neff_s = neff;
%% plot results
for i = 1:nummodes
    neff = neff_(i);
    Ex = (V(1:Nx*Ny,i));
    Ey = (V(Nx*Ny+1:2*Nx*Ny,i));
    Hz = 1i/k0/Z0*(-Uy*Ex + Ux*Ey);
    Hy = erx*Ex/neff/Z0 - 1i*lambdac/2/pi/neff*Uy.'*Hz;
    Hx = -ery*Ey/neff/Z0 - 1i*lambdac/2/pi/neff*Ux.'*Hz;
    Ez = -1i/k0*inv(erz)*Z0*(-Ux.'*Hy + Uy.'*Hx);
    maxfield = max(max(abs(Ex)),max(abs(Ey)))*0.5;
    Exmat = double(zeros(Nx,Ny));
    Eymat = double(zeros(Nx,Ny));
    Ezmat = double(zeros(Nx,Ny));
    Hxmat = double(zeros(Nx,Ny));
    Hymat = double(zeros(Nx,Ny));
    Hzmat = double(zeros(Nx,Ny));
    for y = 1:Ny
        Exmat(:,y) = V(((y-1)*Nx+1:(y)*Nx),i);
        Eymat(:,y) = V(((y-1)*Nx+1:(y)*Nx)+Nx*Ny,i);
        Ezmat(:,y) = imag(Ez((y-1)*Nx+1:y*Nx))/maxfield;
        Hxmat(:,y) = real(Hx((y-1)*Nx+1:y*Nx))/maxfield;
        Hymat(:,y) = real(Hy((y-1)*Nx+1:y*Nx))/maxfield;
        Hzmat(:,y) = imag(Hz((y-1)*Nx+1:y*Nx))/maxfield;
    end
    fighandle = figure(i);
    set(fighandle,'WindowStyle','docked');
    clf
    subplot(221)
    h = surface(x_,y_,ermat.');
    caxis([1 max(max(ermat))])
    set(h,'linestyle','none')
    set(gca,'fontsize',12)
    title('\fontsize{14}Permittivity')
    xlabel('\fontsize{12}x (\mum)')
    ylabel('\fontsize{12}y (\mum)')
    axis([min(x_) max(x_) min(y_) max(y_)])
    view(0,90)
    axis equal
    colormap jet
    colorbar
    subplot(223)
    h = surface(x_,y_,real(Exmat).');
    caxis([-max(max(abs(Exmat))) max(max(abs(Exmat)))])
    set(h,'linestyle','none')
    set(gca,'fontsize',12)
    title(['\fontsize{14}E_x, n_{eff} = ' num2str(neff)])
    xlabel('\fontsize{12}x (\mum)')
    ylabel('\fontsize{12}y (\mum)')
    axis([min(x_) max(x_) min(y_) max(y_)])
    view(0,90)
    axis equal
    subplot(224)
    h = surface(x_,y_,real(Eymat).');
    caxis([-max(max(abs(Eymat))) max(max(abs(Eymat)))])
    set(h,'linestyle','none')
    set(gca,'fontsize',12)
    title(['\fontsize{14}E_y, n_{eff} = ' num2str(neff)])
    xlabel('\fontsize{12}x (\mum)')
    ylabel('\fontsize{12}y (\mum)')
    axis([min(x_) max(x_) min(y_) max(y_)])
    view(0,90)
    axis equal
end