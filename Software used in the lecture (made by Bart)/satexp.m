function [f] = satexp(t,state)

%
% Bullet dynamics is two dimensions
%
xp = state(1);
yp = state(2);
zp = state(3);
xv = state(4);
yv = state(5);
zv = state(6);

GMe = 3.9860044e14;
 
%v = sqrt(xv*xv+yv*yv);
d = sqrt(xp*xp+yp*yp+zp*zp);
v = sqrt(xv*xv+yv*yv+zv*zv);
ad = -GMe./(d*d);

evx = xv/v;
evy = yv/v;
evz = zv/v;

% thrust force
N = 90e-3; % mN PPS-1350 Hall effect

% update forces on model 
xa = ad.* xp./d + N*evx;
ya = ad.* yp./d + N*evy;
za = ad.* zp./d + N*evz;

f = [state(4) state(5) state(6) xa ya za]';