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

%options of the integrator
options = odeset('RelTol', 1e-12);
T = 60000000; %s

[t,y] = ode113(@SunEarthAcc, [0 T], [rE 0 0 0 vE 0 rE+rS 0 0 0 vE+vS 0],options);


figure
hold on
plot3(y(:,1),y(:,2),y(:,3),'b')
plot3(y(:,7),y(:,8),y(:,9),'g')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal
axis vis3d