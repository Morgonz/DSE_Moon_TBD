function [f] = SystemTestAcc(t,state)

%
% Bullet dynamics is two dimensions
%
xp = state(1);
yp = state(2);
zp = state(3);
% xv = state(4);
% yv = state(5);
% zv = state(6);

GMe = 3.9860044e14;
 
%v = sqrt(xv*xv+yv*yv);
d = sqrt(xp*xp+yp*yp+zp*zp);

ad = -GMe./(d*d);

% update forces on model 
xa = ad.* xp./d;
ya = ad.* yp./d;
za = ad.* zp./d;

f = [state(4) state(5) state(6) xa ya za]';