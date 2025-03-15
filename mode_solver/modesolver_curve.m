%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finite difference mode solver
% Chris Doerr 2010
% Bell Labs, Alcatel-Lucent
% all length units in um
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%% adjustable variables
lambdac = 1.55;     % vacuum center wavelength
xlim = 2.0;         % x extent is 2*xlim
ylim = 1.0;         % y extent is 2*ylim
Dx = 0.01;          % x grid size
Dy = 0.01;          % y grid size
Nx = round(xlim*2/Dx);
Ny = round(ylim*2/Dy);
k0 = 2*pi/lambdac;
%% waveguide parameters
R = 3;              % bend radius (RT = R)
wg.w = 1.0;         % waveguide width
wg.h = 0.22;        % waveguide height
ncore = 3.48;       % index of core
nclad = 1.45;       % index of cladding
neff_guess = 2.9;   % initial guess for effective index
%% fill in er
ermat = ones(Nx,Ny)*nclad^2;
%---cladding
for x = -xlim+Dx:Dx:xlim
    ermat(round((x+xlim)/Dx),:) = nclad^2*exp(x/R);
end
%---core
ylim1 = round(-(wg.h/2)/Dy + Ny/2);
ylim2 = round((wg.h/2)/Dy + Ny/2 - 1);
for x = R*log((R-wg.w/2)/R):Dx:R*log((R+wg.w/2)/R)
    ermat(round((x+xlim)/Dx),ylim1:ylim2) = ncore^2*exp(x/R);
end
%% calculate erx, ery, erz
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
display('Setting up matrices...')
temp = ones(Nx*Ny,1);
Ux = spdiags([-temp temp],[0 1],Nx*Ny,Nx*Ny)/Dx;
Uy = spdiags([-temp temp],[0 Nx],Nx*Ny,Nx*Ny)/Dy;
Uy(end,end) = 1/Dy;
P1 = [Ux/erz*Uy.' k0^2*speye(Nx*Ny)-Ux/erz*Ux.'; -k0^2*speye(Nx*Ny)+Uy/erz*Uy.' -Uy/erz*Ux.'];
P2 = [-Ux.'*Uy -k0^2*ery+Ux.'*Ux; k0^2*erx-Uy.'*Uy Uy.'*Ux];
P = P1*P2/k0^4;
%% find eigenvalues
display('Calculating eigenvalues...')
options.disp = 0;
options.tol = 1e-9;
[V,D] = eigs(P,4,neff_guess^2,options);
neff = sqrt(diag(D));
display(neff)
%% plot waveguide structure and modes
clf
colormap jet
modenum = 1;   % mode number to plot. You can change it and rerun this cell
for y = 1:Ny;
    Exmat(:,y) = V(((y-1)*Nx+1:(y)*Nx),modenum);
    Eymat(:,y) = V(((y-1)*Nx+1:(y)*Nx)+Nx*Ny,modenum);
end
x = ((1:Nx) - Nx/2)*Dx;
y = ((1:Ny) - Ny/2)*Dy;
subplot(221)
h = surface(x,y,real(ermat).');
caxis([-max(max(real(ermat))) max(max(real(ermat)))])
set(h,'linestyle','none')
title('\fontsize{14}Permittivity')
xlabel('\fontsize{12}x (\mum)')
ylabel('\fontsize{12}y (\mum)')
view(0,90)
subplot(223)
h = surface(x,y,real(Exmat).');
caxis([-max(max(abs(Exmat))) max(max(abs(Exmat)))])
set(h,'linestyle','none')
title(['\fontsize{14}E_x, n_{eff} = ' num2str(neff(modenum))])
xlabel('\fontsize{12}x (\mum)')
ylabel('\fontsize{12}y (\mum)')
view(0,90)
subplot(224)
h = surface(x,y,real(Eymat).');
caxis([-max(max(abs(Eymat))) max(max(abs(Eymat)))])
set(h,'linestyle','none')
title(['\fontsize{14}E_y, n_{eff} = ' num2str(neff(modenum))])
xlabel('\fontsize{12}x (\mum)')
ylabel('\fontsize{12}y (\mum)')
view(0,90)
shg