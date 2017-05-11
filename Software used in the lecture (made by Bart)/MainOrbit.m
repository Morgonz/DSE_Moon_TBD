% This program is a test program to see the applications of the animation
% of a cannon ball

clear all;
close all;
clc;

% Initial conditions
GMe = 3.9860044e14;
Re =  6371000;
h = 500000;

% circular orrbit: initial conditions
V0 = sqrt(GMe./(Re+h)); % m/s 
T = 2*pi*sqrt((Re+h)^3/GMe); % sec

% Transfer orbit: initial conditions
r2 = 385000600;
r1 = Re+h;

aT = (r1+r2)/2;

DV1 = sqrt(GMe/r1)*(sqrt((2*r2)/(r1+r2))-1);
DV2 = sqrt(GMe/r2)*(1-sqrt((2*r1)/(r1+r2)));
T2 = pi*sqrt((aT)^3/GMe);

% options of the integrator
options = odeset('RelTol',1e-12);

% The time integration of the sat
[t,y] = ode45(@satacc,[0 T],[0 Re+h 0 V0 0 0],options);
[t2,y2] = ode45(@satacc,[0 T2],[0 Re+h 0 V0+DV1 0 0],options);

% earth
[X,Y,Z] = sphere(40);

figure
plot3(y(:,1),y(:,2),y(:,3),'LineWidth',4)
hold on
surf(X*Re, Y*Re, Z*Re) % where (a,b,c) is center of the sphere
plot3(y2(:,1),y2(:,2),y2(:,3),'r')
hold off
axis equal
