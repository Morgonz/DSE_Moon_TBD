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
theta = -180; %deg angle the satellite is at the start of the simulation with respect to the positive x-axis

%orbit of the Moon
rM = 385000600; %m
vM = sqrt(Mu_Earth/rM); %m/s speed of rotation of the Moon

%orbit at the Moon
h3 = 9000000; %m
r3 = R_Moon + h3; %m
v3 = sqrt(Mu_Moon/r3); %m/s


%options of the integrator
options1 = odeset('RelTol', 1e-9);
options2 = odeset('RelTol', 1e-9, 'Events', @CrossMoonOrbit);
T = 8640000; %s time of simulation

%[t,y] = ode45(@StartingOrbitAcc,[0 10*T],[-201642375	237221957.5	0	-606.4649404	-1043.191023	0	-247843716.9	284494001.5	0	-776.6183011	-689.085314	0],options1); %Starting orbit of the satellite
[t,y,te] = ode45(@LowThrustAcc,[0 4*T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V0*sin(theta*(pi/180)) V0*cos(theta*(pi/180))  0 rM 0 0 0 vM 0],options2); %Transfer orbit
[t3,y3] = ode45(@StartingOrbitAcc,[0 T/5],[y(end,1) y(end,2) y(end,3) y(end,4) y(end,5) y(end,6) y(end,7) y(end,8) y(end,9) y(end,10) y(end,11) y(end,12)],options1); %Starting orbit of the satellite
%defining the sphere of the Earth
[X,Y,Z] = sphere(40);

%plotting the figure
figure
plot3(y(:,1),y(:,2),y(:,3),'b') %Low-Earth orbit of the satellite
hold on
plot3(y(end,1),y(end,2),y(end,3),'rh') %mark last location of the satellite
plot3(y(:,7),y(:,8),y(:,9),'r') %orbit of the Moon 
plot3(y(1,7),y(1,8),y(1,9),'r*') %mark starting location of the Moon
plot3(y3(end,7),y3(end,8),y3(end,9),'ro') %mark last location of the Moon
plot3(y3(:,1),y3(:,2),y3(:,3),'c') %orbit during coasting
plot3(y3(:,7),y3(:,8),y3(:,9),'r')%Orbit of the Moon during coasting
%surf(X*R_Earth, Y*R_Earth, Z*R_Earth) % where (a,b,c) is center of the sphere
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal


%Define Event Function
function [value,isterminal,direction] = CrossMoonOrbit(t,y)
value = sqrt((y(1))^2 + (y(2))^2 + (y(3))^2)-385000600+66100000;
isterminal = 1;
direction = 1;
end


