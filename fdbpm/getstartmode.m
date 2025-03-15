% function that calculates the waveguide mode using the effective index method
% matches the beam-propagation program

function [Ey] = getstartmode(n,nc,w,lambda,x)

k = 2*pi/lambda;

V = k*w*sqrt(n^2 - nc^2);

bnorm = 0; nu = 0;
while V > (nu*pi+2*atan(sqrt(bnorm/(1-bnorm))))/sqrt(1-bnorm)
   bnorm = bnorm + 0.0005;		
   if bnorm>1
      bnorm = 1;
      break
      end
end

b = k*sqrt(nc^2 + bnorm*(n^2 - nc^2));

N = 1;

h = w;

g = sqrt(-(nc*k)^2 + b^2);

a = sqrt((n*k)^2 - b^2);

A = zeros(3,3);
A(1,1) = 1;
A(2,2) = a;

v = zeros(3,1);
v(1) = 1;
v(2) = g;

A(3,1) = cos(a*h);
A(3,2)   = sin(a*h);
A(3,3) = -exp(-g*h);

u = inv(A)*v;

N = max(size(x));
for cnt = 1:N
   if x(cnt) <= -h/2
      Ey(cnt) = exp(g*(x(cnt) + h/2));
   elseif x(cnt) >= h/2
      Ey(cnt) = u(3)*exp(-g*(x(cnt) + h/2));
   else
      Ey(cnt) = u(1)*cos(a*(x(cnt) + h/2)) + u(2)*sin(a*(x(cnt) + h/2));
   end
end

