function [f] = sat3BP(t,state)

%
% Bullet dynamics is two dimensions
%
xp = state(1);
yp = state(2);
zp = state(3);
% xv = state(4);
% yv = state(5);
% zv = state(6);

xM = state(7);
yM = state(8);
zM = state(9);

GMe = 3.9860044e14;
GMm = 4.902801e12;
 
%v = sqrt(xv*xv+yv*yv);
d = sqrt(xp*xp+yp*yp+zp*zp);
dm = sqrt((xM-xp)^2+(yM-yp)^2+(zM-zp)^2);
dM = sqrt(xM*xM+yM*yM+zM*zM);

ad = -GMe./(d*d*d);
adM = -GMm./(dm*dm*dm);
aM = -(GMe+GMm)./(dM*dM*dM);

% update forces on model 
xa = ad.* xp + adM*(xp-xM);
ya = ad.* yp + adM*(yp-yM);
za = ad.* zp + adM*(zp-zM);

xaM = aM.* xM;
yaM = aM.* yM;
zaM = aM.* zM;

f = [state(4) state(5) state(6) xa ya za state(10) state(11) state(12) xaM yaM zaM]';