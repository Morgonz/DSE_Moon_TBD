clear all;
close all;
clc;

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Moon = 4.902801e12; %m^3 s^-2
R_Earth = 6371000; %m
R_Moon = 1737000; %m

%Circular Orbit initial conditions
h = 500000; %m
V0 = sqrt(Mu_Earth/(R_Earth + h)); %m/s
r1 = R_Earth + h; %m
T1 = 2*pi*sqrt((r1^3)/Mu_Earth);

t_days = (T1*1000)/(3600*24)
 

options1 = odeset('RelTol', 2.2205e-14);
[t,y] = ode113(@SystemTestAcc,[0 1000*T1],[r1 0  0 0 V0  0 ],options1); %initial orbit

figure
hold on
plot3(y(:,1),y(:,2),y(:,3),'b') %initial orbit
plot3(y(end,1),y(end,2),y(end,3),'*k') %initial orbit
axis equal
xlabel('x [m]')
ylabel('y [m]')
title('Orbit Earth')
hold off

R_error = sqrt(y(end,1)^2+y(end,2)^2)
error = (R_error-r1)/r1
error_percentage = 100*error
