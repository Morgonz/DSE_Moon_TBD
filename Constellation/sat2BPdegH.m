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

% Moon inertial to moon fixed
%add rotation to satellites
ome = 2*pi/(29.5*3600*24); %rad/s - https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
%% Time conversion from seconds to juliandate
t_rot = t;
OME  = -ome*t_rot;
Tz_E = rotationmatrices.Tz(OME);
Ty_E = rotationmatrices.Ty(deg2rad(6.69)); % from https://ssd.jpl.nasa.gov/dat/lunar_cmd_2005_jpl_d32296.pdf

%%      
% Earth influence
GMe = 3.9860044e14;
xe = 385000600;
ye = 0;
ze = 0;

% Moon SH
SatM_coords  = Tz_E*[xp yp zp]'; 
Earth_coords = Ty_E*Tz_E*[xe ye ze]';
model        = 'LP165P';
degree       = 50;

de  = sqrt((Earth_coords(1)-xp)^2+(Earth_coords(2)-yp)^2+(Earth_coords(3)-zp)^2);
dE  = sqrt((Earth_coords(1))^2+(Earth_coords(2))^2+(Earth_coords(3))^2);
ade = -GMe./(de*de*de); %add moon gravity loc
adE = -GMe./(dE*dE*dE); 

[gx, gy, gz] = gravitysphericalharmonic(SatM_coords',model,degree);
G_arr = Tz_E'*[gx, gy, gz]';

% minus due to vector from Moon (maybe it doens't work since this actually
% implies there is also a VM (more momentum less the influence of forces)
gEx = -ade*(Earth_coords(1)-xp)-(-adE*Earth_coords(1));
gEy = -ade*(Earth_coords(2)-yp)-(-adE*Earth_coords(2));
gEz = -ade*(Earth_coords(3)-zp)-(-adE*Earth_coords(3));
% update forces on model 
% Attract due Moon and Earth
xa = G_arr(1) + gEx;
ya = G_arr(2) + gEy;
za = G_arr(3) + gEz;

f = [state(4) state(5) state(6) xa ya za]';