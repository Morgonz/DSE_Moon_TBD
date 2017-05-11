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
T = 1* 24 * 3600; % sec

% options of the integrator
options = odeset('RelTol',1e-6);

% The time integration of the sat
[t,y] = ode45(@satacc,[0 T],[0 Re+h 0 V0 0 0],options);
[t2,y2] = ode45(@satexp,[0 T],[0 Re+h 0 V0 0 0],options);

figure
plot3(y(:,1),y(:,2),y(:,3))
hold on
plot3(y2(:,1),y2(:,2),y2(:,3),'g')
hold off
axis equal
