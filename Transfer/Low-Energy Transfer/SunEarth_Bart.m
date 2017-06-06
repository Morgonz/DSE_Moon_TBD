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

%Satellite initial conditions
h = 500000; %m
vS = sqrt(Mu_Earth/(R_Earth + h)); %m/s
rS = R_Earth + h; %m

%Sun-Earth L2
rL = rE*(Mu_Earth/(3*Mu_Sun))^(1/3);

%options of the integrator
options = odeset('RelTol', 1e-18);
T = 10000000; %s

[t,y] = ode113(@SunEarthAcc, [0 T], [rE 0 0 0 vE 0 rE+rL-200000000 0 0 0 vE+vS/10 0],options);

yrot = RotatingFrameSunEarth(y);

figure
hold on
plot3(yrot(:,7),yrot(:,8),yrot(:,9))
plot3(rE+rL,0,0,'c*')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal
axis vis3d