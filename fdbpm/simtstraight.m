function [nz,nx] = simtstraight(lstart,lstop,a,wstart,wstop,z0,x0)

global Dz
global w
global d

hstart = round((z0 + lstart*cos(a))/Dz) + 1;
hstop = round((z0 + lstop*cos(a))/Dz) + 1;

for h = hstart:hstop
    z = h*Dz;
    w(h) = (wstop + (wstart - wstop)/(hstop - hstart)*(hstop - h))/cos(a);
    d(h) = (z - z0)*tan(a) + x0;
end

nz = z0 + lstop*cos(a);
nx = x0 + lstop*sin(a);

