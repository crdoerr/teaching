function [ermat] = addwgpiece(x1,y1,x2,y2,n,ermat,Dx,Dy)
% adds piece to permittivity distribution
% the points have their origin at the center of the distribution
% Chris Doerr
[Nx,Ny] = size(ermat);
xlim1 = round(x1/Dx + Nx/2);
xlim2 = round(x2/Dx + Nx/2 - 1);
ylim1 = round(y1/Dy + Ny/2);
ylim2 = round(y2/Dy + Ny/2 - 1);
%% clip to limits
xlim1 = bound(xlim1,1,Nx);
xlim2 = bound(xlim2,xlim1,Nx);
ylim1 = bound(ylim1,1,Ny);
ylim2 = bound(ylim2,ylim1,Ny);
%% load array
for y = ylim1:ylim2
    ermat(xlim1:xlim2,y) = n^2;
end

function y = bound(x,bl,bu)
  % return bounded value clipped between bl and bu
  y=min(max(x,bl),bu);