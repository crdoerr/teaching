function [nz,nx] = simcurve(type,astart,astop,R,Cw,Co,z0,x0)

global Dz
global w
global d

if type=='L'
    R = -R;
    Co = -Co;
    Cw = -Cw;
end

zorigin = R*sin(astart) + z0;
xorigin = -R*cos(astart) + x0;

hstart = round(z0/Dz) + 1;
hstop = round((z0 + R*(sin(astart) - sin(astop)))/Dz) + 1;

for h = hstart:hstop
    z = h*Dz;
    a = asin((zorigin - z)/R);
    w(h) = abs(Cw/cos(a));
    if type=='R'
        d(h) = xorigin + sqrt((R-Co)^2 - (z - zorigin)^2);
    else
        d(h) = xorigin - sqrt((R-Co)^2 - (z - zorigin)^2);
    end
end

nz = z0 + R*sin(astart) - R*sin(astop);
nx = x0 - R*cos(astart) + R*cos(astop);

