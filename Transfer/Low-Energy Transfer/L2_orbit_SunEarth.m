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

% SUN
[X,Y,Z] = sphere(40);

%options of the integrator
options1 = odeset('RelTol', 1e-18);

%time option integrator
T_days = 500;
T = 60*60*24*T_days; %s
t_kick_days = 10
t_kick_sec = 60*60*24*t_kick_days

[t1,y1] = ode113(@SunEarthAcc, [0 T], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options1);
yrot1 = RotatingFrameSunEarth(y1);
%{
% obtain end velocity
V_abs = sqrt(y1(end,10)^2+y1(end,11)^2)
kick = 10 %m/s
dv_x = (y1(end,10)/V_abs)*kick
dv_y = (y1(end,11)/V_abs)*kick

[t2,y2] = ode113(@SunEarthAcc, [0 T-t_kick_sec],[y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+dv_x y1(end,11)+dv_y y1(end,12)],options1);
yrot2 = RotatingFrameSunEarth(y2);
%}
hold on
%start plotting
plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9),'DisplayName','Halo Orbit');
%plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9),'DisplayName','Halo Orbit part 2');
L2 = plot3(rE+rL,0,0,'k*','DisplayName','Earth-Sun L2');
Earth = plot3(rE,0,0,'bo','DisplayName','Earth');
%surf(X*R_Sun,Y*R_Sun, Z*R_Sun,'DisplayName','Sun position')
hold off
title('L2 Orbit SunEarth')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
%legend([L2 Earth],{'Earth-Sun L2','Earth'})
legend('show')
axis equal
axis vis3d