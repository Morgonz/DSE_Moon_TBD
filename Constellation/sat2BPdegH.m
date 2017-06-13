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
ETB = 0;
STB = 0;
MIR = 1;

t_rot = t;
%% Earth Third Body Perturbations
if ETB==1
    ome_E = 2*pi/(27.32*3600*24); %rad/s - https://ssd.jpl.nasa.gov/dat/lunar_cmd_2005_jpl_d32296.pdf
    OME_E  = -ome_E*t_rot;
    
    Tz_E = rotationmatrices.Tz(OME_E);
    Ty_E = rotationmatrices.Ty(deg2rad(6.69)); % from https://ssd.jpl.nasa.gov/dat/lunar_cmd_2005_jpl_d32296.pdf

    GMe = 3.9860044e14;
    xe = 365000600; ye = 0; ze = 0; %
    
    E_coords     = Ty_E*Tz_E*[xe ye ze]';
    
    de  = sqrt((E_coords(1)-xp)^2+(E_coords(2)-yp)^2+(E_coords(3)-zp)^2);
    dE  = sqrt((E_coords(1))^2+(E_coords(2))^2+(E_coords(3))^2);
    ade = -GMe./(de*de*de); 
    adE = -GMe./(dE*dE*dE); 
    
    gEx = -ade*(E_coords(1)-xp)-(-adE*E_coords(1));
    gEy = -ade*(E_coords(2)-yp)-(-adE*E_coords(2));
    gEz = -ade*(E_coords(3)-zp)-(-adE*E_coords(3));
else
    gEx = 0; 
    gEy = 0;
    gEz = 0;
end

if STB==1
    ome_S = 2*pi/(365*3600*24); %rad/s - https://ssd.jpl.nasa.gov/dat/lunar_cmd_2005_jpl_d32296.pdf
    OME_S  = -ome_S*t_rot;

    Tz_S = rotationmatrices.Tz(OME_S);
    
    GMs = 1.32712e20;
    xs = 149.60e9; ys = 0; zs = 0; 
    
    S_coords     = Tz_S*[xs ys zs]';
    
    ds  = sqrt((S_coords(1)-xp)^2+(S_coords(2)-yp)^2+(S_coords(3)-zp)^2);
    dS  = sqrt((S_coords(1))^2+(S_coords(2))^2+(S_coords(3))^2);
    ads = -GMs./(ds*ds*ds); 
    adS = -GMs./(dS*dS*dS);
    
    gSx = -ads*(S_coords(1)-xp)-(-adS*S_coords(1));
    gSy = -ads*(S_coords(2)-yp)-(-adS*S_coords(2));
    gSz = -ads*(S_coords(3)-zp)-(-adS*S_coords(3));
else 
    gSx = 0;
    gSy = 0;
    gSz = 0;
end

%% Moon influence  
if MIR==1
    ome_E = 2*pi/(27.32*3600*24); %rad/s - https://ssd.jpl.nasa.gov/dat/lunar_cmd_2005_jpl_d32296.pdf
    OME_E  = -ome_E*t_rot;

    Tz_E = rotationmatrices.Tz(OME_E);
    
    model        = 'LP165P';
    degree       = 50;
    SatM_coords  = Tz_E*[xp yp zp]'; 

    [gx, gy, gz] = gravitysphericalharmonic(SatM_coords',model,degree);
    
    G_arr = Tz_E'*[gx, gy, gz]';
else
    GMm = 4.904e12;
    
    dm  = sqrt((xp)^2+(yp)^2+(zp)^2);
    adm = -GMm./(dm*dm*dm);
    G_arr =[ adm*(xp) adm*(yp) adm*(zp)];
end
    
% update forces on model 
% Attract due Moon and Earth
xa = G_arr(1) + gEx + gSx;
ya = G_arr(2) + gEy + gSy;
za = G_arr(3) + gEz + gSz;

f = [state(4) state(5) state(6) xa ya za]';
