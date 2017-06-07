clear all;
close all;
clc;

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Sun = 1.327124E+20; %m^3 s^-2 from wikipedia
R_Earth = 6371000; %m
R_Sun = 695508000; %m https://solarsystem.nasa.gov/planets/sun/facts

%Earth initial conditions
rE = 149.60E+9; %m radius of the Earth https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
vE = sqrt((Mu_Earth+Mu_Sun)/rE); %rotational speed of the Earth


%Sun-Earth L2
rL = rE*(Mu_Earth/(3*Mu_Sun))^(1/3);

%options of the integrator
options = odeset('RelTol', 2.22045e-14);
days = 800;
T = 60*60*24*days; %s

[t,y] = ode113(@SunEarthAcc, [0 T], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options);
options1 = odeset('RelTol', 2.22045e-14, 'Events', @LeaveHalo);

%{
yrot = RotatingFrameSunEarth(y);
figure
hold all
title(['Used dv is:'])
plot3(yrot(:,7),yrot(:,8),yrot(:,9))
plot3(rE+rL,0,0,'c*')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal
axis vis3d
%}
%
for dv = linspace(-1000,1000,3)
%disp(dv)

[t,y,te] = ode113(@SunEarthAcc, [0 T], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426+dv 1.260005700065],options1);
if te>2500
    disp(['ending time for dv of ' int2str(dv) ' is:'])
    disp(te)   
else
    yrot = RotatingFrameSunEarth(y);
    figure
    hold all
    title(['Used dv is:' int2str(dv)])
    plot3(yrot(:,7),yrot(:,8),yrot(:,9))
    plot3(rE+rL,0,0,'c*')
    hold off
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    axis equal
    axis vis3d
end
end
%}
function [value,isterminal,direction] = LeaveHalo(t,y)
value = 3*10^8 - sqrt((y(7)-1.509765473987467e+11)^2 + (y(8))^2) ;
isterminal = 1;
direction = 0;
end