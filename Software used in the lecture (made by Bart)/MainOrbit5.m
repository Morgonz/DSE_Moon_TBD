% This program is a test program to see the applications of the animation
% of a cannon ball

clear all;
close all;
clc;

% Initial conditions
GMe = 3.9860044e14;
GMm = 4.902801e12;
Re =  6371000;
h = 500000;
r2 = 385000600;

% circular orrbit: initial conditions
% circular orrbit: initial conditions
V0 = sqrt(GMe./(Re+h)); % m/s 
%T = 2*pi*sqrt((Re+h)^3/GMe); % sec
T = 20* 24 * 3600; % sec

% Transfer orbit: initial conditions
r2 = 385000600;
R1 = r2*(GMm/(3*GMe))^(1/3);
VM = sqrt(GMe./(r2)); % m/s 
r1 = Re+h;

DV1 = sqrt(GMe/r1)*(sqrt((2*(r2-R1))/(r1+(r2-R1)))-1);

aT = (r1+r2-R1)/2;
T2 = 40*pi*sqrt((aT)^3/GMe);


% options of the integrator
options = odeset('RelTol',1e-12);

% The time integration of the sat
[t,y] = ode45(@sat3BP,[0 T2],[0 Re+h 0 V0+DV1 0 0 0 r2 0 VM 0 0],options);

figure
subplot(1,2,1)
plot3(y(:,1),y(:,2),y(:,3))
hold on
plot3(y(:,7),y(:,8),y(:,9),'r')
plot3(0,0,0,'ob')
hold off
axis equal

% rotational

theta = atan2(y(:,8),y(:,7));

xr = cos(theta).*y(:,1)  + sin(theta).*y(:,2);
yr = -sin(theta).*y(:,1) + cos(theta).*y(:,2);
zr = y(:,3);

xrM = cos(theta).*y(:,7)  + sin(theta).*y(:,8);
yrM = -sin(theta).*y(:,7) + cos(theta).*y(:,8);
zrM = y(:,9);

subplot(1,2,2)
plot3(xr./1e3,yr./1e3,zr./1e3)
hold on
plot3(xr(1)./1e3,yr(1)./1e3,zr(1)./1e3,'*')
plot3(xrM(1)./1e3,yrM(1)./1e3,zrM(1)./1e3,'or')
plot3(0,0,0,'ob')
hold off
axis equal