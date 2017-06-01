clear all;
close all;
clc;

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Moon = 4.902801e12; %m^3 s^-2
R_Earth = 6371000; %m mean radius of the Earth
R_Moon = 1737400; %m mean Moon radius

%Circular Orbit initial conditions
h = 500000; %m
V0 = sqrt(Mu_Earth/(R_Earth + h)); %m/s
r1 = R_Earth + h; %m

%orbit of the Moon
rM = 385000600; %m
vM = sqrt((Mu_Earth+Mu_Moon)/rM); %m/s speed of rotation of the Moon

%Transfer orbit
DV1 = sqrt(Mu_Earth/r1)*(sqrt((2*rM)/(r1+rM))-1);
theta=-124;
V1 = V0+DV1;


%options of the integrator
options = odeset('RelTol', 1e-12);
options1 = odeset('RelTol', 1e-12, 'Events', @CrossMoonOrbit);
T=10000000;

%integrator
[t,y] = ode45(@HohmannAcc,[0 T],[0 r1 0 -V0 0 0 rM 0 0 0 vM 0],options); %initial orbit
[t2,y2,te] = ode45(@HohmannAcc,[0 T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V1*sin(theta*(pi/180)) V1*cos(theta*(pi/180)) 0 rM 0 0 0 vM 0],options1); %Transfer Orbit

[X,Y,Z] = sphere(20);
figure
hold on
plot3(y2(:,1),y2(:,2),y2(:,3));
plot3(y2(:,7),y2(:,8),y2(:,9));
surf(X*R_Earth,Y*R_Earth,Z*R_Earth);
surf(X*R_Moon+y2(end,7),Y*R_Moon+y2(end,8),Z*R_Moon+y2(end,9));
hold off
axis equal
axis vis3d
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Inertial Frame')

yrot = RotatingFrame(y2);

figure
hold on
plot3(yrot(:,1),yrot(:,2),yrot(:,3));
plot3(yrot(:,7),yrot(:,8),yrot(:,9));
surf(X*R_Earth,Y*R_Earth,Z*R_Earth);
surf(X*R_Moon+yrot(end,7),Y*R_Moon+yrot(end,8),Z*R_Moon+yrot(end,9));
hold off
axis equal
axis vis3d
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Rotating Frame')

%Define Event Function, target orbit at 1000 km above Moon surface
function [value,isterminal,direction] = CrossMoonOrbit(t2,y2)
value = sqrt((y2(1)-y2(7))^2 + (y2(2)-y2(8))^2 + (y2(3)-y2(9))^2)-2436000;
isterminal = 1;
direction = 0;
end
