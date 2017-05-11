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
%T = 2*pi*sqrt((Re+h)^3/GMe); % sec
T = 20* 24 * 3600; % sec

% Transfer orbit: initial conditions
r2 = 385000600;
VM = sqrt(GMe./(r2)); % m/s 
r1 = Re+h;

aT = (r1+r2)/2;

DV1 = sqrt(GMe/r1)*(sqrt((2*r2)/(r1+r2))-1);
DV2 = sqrt(GMe/r2)*(1-sqrt((2*r1)/(r1+r2)));
T2 = pi*sqrt((aT)^3/GMe);

% options of the integrator
options = odeset('RelTol',1e-6);

% The time integration of the sat
[t,y] = ode45(@sat3BP,[0 T],[0 Re+h 0 V0 0 0 0 -r2 0 -VM 0 0],options);
[t2,y2] = ode45(@sat3BP,[0 T2],[0 Re+h 0 V0+DV1 0 0 0 -r2 0 -VM 0 0],options);

figure
plot3(y(:,1),y(:,2),y(:,3))
hold on
plot3(y(:,7),y(:,8),y(:,9),'r')
plot3(y2(:,1),y2(:,2),y2(:,3),'g')
hold off
axis equal
