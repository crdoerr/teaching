%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-D FD-BPM simulator
% Chris Doerr
% Copyright 2011
% Not for commercial use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
global Dz
global w
global d
w = []; d = []; d1 = []; d2 = [];
%% structure parameters
lambdac = 1.546;        % design wavelength
ncore = 2.8;            % core effective index
nclad = 1.45;           % cladding index
SMFw = 1.0;           % width of single-mode waveguide
%% directional coupler parameters
dc.w = 0.4;
dc.sp = 0.68;
dc.l = 5;
dc.R = 25;
dc.a = 8*pi/180;
%% simulation parameters
wstart = dc.w;        % width of starting waveguide
N = 2;                 % step for recording data
lambda = 1.546;         % simulation wavelength
Dz = 0.25;              % step size in propagation direction
numpts = 2^12;          % number of points in transverse direction
xlimit = 8;            % window extent (xlimit = half window size)
n = ncore;
Dx = 2*xlimit/numpts;
Dkx = 2*pi/Dx/numpts;
x = -xlimit:Dx:xlimit-Dx;
kx = -pi/Dx:Dkx:pi/Dx-Dkx;
k = n*2*pi/lambda;
kz = sqrt(k.^2 - kx.^2);
%% put structure into w and d
[nz,nx] = simstraight(0,10,0,dc.w,0,0);
[nz,nx] = simcurve('R',0,-dc.a,dc.R,dc.w,0,nz,nx);
[nz,nx] = simstraight(0,10,-dc.a,dc.w,nz,nx);
[nz,nx] = simcurve('L',-dc.a,0,dc.R,dc.w,0,nz,nx);
dcoup = d(end);
[nz,nx] = simstraight(0,dc.l,0,dc.w,nz,nx);
[nz,nx] = simcurve('L',0,dc.a,dc.R,dc.w,0,nz,nx);
[nz,nx] = simstraight(0,10,dc.a,dc.w,nz,nx);
[nz,nx] = simcurve('R',dc.a,0,dc.R,dc.w,0,nz,nx);
[nz,nx] = simstraight(0,10,0,dc.w,nz,nx);
zcntmax = max(size(w));
z = ((1:zcntmax) - 1)*Dz;
%% find starting, uncoupled waveguide mode
u = getstartmode(ncore,nclad,wstart/4,lambda,x);
u = u/max(abs(u));
%% do FD BPM
pwrin = sum(abs(u).^2);
kslowvar = 2*pi*ncore/lambda;
Z = ones(round(length(w)/10),length(x));
Y = Z;
alpha = zeros(1,numpts);
beta = zeros(1,numpts);
unew = zeros(1,numpts);
for cnt = 1:length(w)
    k = 2*pi*nclad/lambda*ones(size(u));
    lim1 = round((d(cnt) - w(cnt)/2)/Dx) + numpts/2 + 1;
    lim2 = round((d(cnt) + w(cnt)/2)/Dx) + numpts/2 + 1;
    k(lim1:lim2) = 2*pi*ncore/lambda;
    lim1 = round((-d(cnt) - w(cnt)/2 + 2*dcoup - dc.sp)/Dx) + numpts/2 + 1;
    lim2 = round((-d(cnt) + w(cnt)/2 + 2*dcoup - dc.sp)/Dx) + numpts/2 + 1;
    k(lim1:lim2) = 2*pi*ncore/lambda;
    s = 2 - Dx^2*(k.^2 - kslowvar^2) + 1i*4*kslowvar*Dx^2/Dz;
    q = -2 + Dx^2*(k.^2 - kslowvar^2) + 1i*4*kslowvar*Dx^2/Dz;
    b = s;
    dd = [u(2:end) 0] + q.*u + [0 u(1:end-1)];
    alpha(1) = 1/b(1);
    beta(1) = dd(1)/b(1);
    for i = 2:numpts-1
        alpha(i) = 1/(b(i) - alpha(i-1));
        beta(i) = (dd(i) + beta(i-1))/(b(i) - alpha(i-1));
    end
    unew(numpts) = 0;
    unew(numpts-1) = beta(numpts-1);
    for i = numpts-2:-1:2
        unew(i) = alpha(i)*unew(i+1) + beta(i);
    end
    unew(1) = 0;
    u = unew;
    if cnt/N==round(cnt/N)
        Z(cnt/N,:) = abs(u).^2;
        Y(cnt/N,:) = k - kslowvar;
    end
end
pwrout = sum(abs(u).^2);
%% plot results
fprintf('Power is conserved to within %4.3f percent.\n',abs((pwrout-pwrin)/pwrin)*100)
clf
sizeZ = size(Z);
m = 1:sizeZ(1);
h = surface(m*Dz*N,x,(Z.'-Y.'/100));
% caxis([0 0.5])
colormap('jet')  %'colorcube' is also nice
colormap(flipud(colormap))
view(0,90)
set(h,'linestyle','none')
axis([0,max(z),-xlimit,xlimit])
set(gca,'fontsize',12)
xlabel('\fontsize{12}z (\mum)')
ylabel('\fontsize{12}x (\mum)')
shg