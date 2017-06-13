function [f] = sat2BP(t,state)

%
% Bullet dynamics is two dimensions
%
xp = state(1);
yp = state(2);
zp = state(3);
% xv = state(4);
% yv = state(5);
% zv = state(6)
MIR = 0;
if MIR
    Moon_coords = [xp yp zp]; 
    model       = 'LP165P';
    degree      = 0;

    [gx, gy, gz]  = gravitysphericalharmonic(Moon_coords,model,degree);
    G_arr = [gx gy gz];
else
    GMm = 4.904e12;
    dm  = sqrt((xp)^2+(yp)^2+(zp)^2);
    adm = -GMm./(dm*dm*dm);
    G_arr =[ adm*(xp) adm*(yp) adm*(zp)];
end

xa = G_arr(1);
ya = G_arr(2);
za = G_arr(3);

f = [state(4) state(5) state(6) xa ya za]';