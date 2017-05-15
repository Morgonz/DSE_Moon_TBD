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
theta = -19; %deg angle the satellite is at the start of the simulation with respect to the positive x-axis

%orbit of the Moon
rM = 385000600; %m
vM = sqrt(Mu_Earth/rM); %m/s speed of rotation of the Moon

%options of the integrator
options1 = odeset('RelTol', 1e-12);
options2 = odeset('RelTol', 1e-12, 'Events', @CrossMoonOrbit);
T = 8640000; %s time of simulation

[t,y] = ode45(@StartingOrbitAcc,[0 T],[0 r1 0 -V0 0 0 rM 0 0 0 vM 0],options1); %Starting orbit of the satellite
[t2,y2,te] = ode45(@LowThrustAcc,[0,T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V0*sin(theta*(pi/180)) V0*cos(theta*(pi/180))  0 rM 0 0 0 vM 0],options2); %Transfer orbit

%defining the sphere of the Earth
[X,Y,Z] = sphere(40);

%plotting the figure
figure
plot3(y(:,1),y(:,2),y(:,3),'b') %Low-Earth orbit of the satellite
hold on
plot3(y2(:,7),y2(:,8),y2(:,9),'r') %orbit of the Moon 
plot3(y2(1,7),y2(1,8),y2(1,9),'r*') %mark starting location of the Moon
plot3(y2(length(y2),7),y2(length(y2),8),y2(length(y2),9),'ro') %mark last location of the Moon
plot3(y2(:,1),y2(:,2),y2(:,3),'g') %Transfer orbit
%surf(X*R_Earth, Y*R_Earth, Z*R_Earth) % where (a,b,c) is center of the sphere
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal

%Define Event Function
function [value,isterminal,direction] = CrossMoonOrbit(t2,y2)
value = sqrt((y2(1)-y2(7))^2 + (y2(2)-y2(8))^2 + (y2(3)-y2(9))^2)-66100000;
isterminal = 1;
direction = 1;
end


