function [ermat] = addwgpiece(x1,y1,x2,y2,n,ermat,x,y)
% adds piece to permittivity distribution
% the points have their origin at the center of the distribution
% Chris Doerr

ermat(x >= x1 & x <= x2, ...
    y >= y1 & y <= y2) = n^2; % set each coordinate's permittivity to n^2

end