function [f] = sat2BP(t,state)

%
% Bullet dynamics is two dimensions
%
xp = state(1);
yp = state(2);
zp = state(3);
% xv = state(4);
% yv = state(5);
% zv = state(6);

Moon_coords = [xp yp zp]; 
model       = 'LP165P';
degree      = 0;

[gx, gy, gz]  = gravitysphericalharmonic(Moon_coords,model,degree);

% update forces on model 
% attract due earth, then attract due moon
xa = gx;
ya = gy;
za = gz;

f = [state(4) state(5) state(6) xa ya za]';