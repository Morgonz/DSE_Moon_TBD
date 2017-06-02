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

%orbit of the Moon
rM = 385000600; %m
vM = sqrt((Mu_Earth+Mu_Moon)/rM); %m/s speed of rotation of the Moon

%transfer orbit
aT  = (r1+rM)/2; 
DV1 = sqrt(Mu_Earth/r1)*(sqrt((2*rM)/(r1+rM))-1); %First Delta V
T = 864000; %s time of simulation
V1 = V0+DV1; %velocity just after the initial burn
theta = -125.2; %deg angle of initial position of satellite with respect to positive x-axis

%insertion in the orbit
xDV = -1026.13455674416;
yDV = -412.741348298507;
zDV = 906.131403670556;
DV2 = sqrt(xDV^2 + yDV^2 + zDV^2);

%options of the integrator
options = odeset('RelTol', 1e-12);
options1 = odeset('RelTol', 1e-12, 'Events', @CrossMoonOrbit);

[t,y] = ode45(@HohmannAcc,[0 T],[0 r1 0 -V0 0 0 rM 0 0 0 vM 0],options); %initial orbit
[t2,y2,te] = ode45(@HohmannAcc,[0 T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V1*sin(theta*(pi/180)) V1*cos(theta*(pi/180)) 0 rM 0 0 0 vM 0],options1); %Transfer Orbit
[t3,y3] = ode45(@HohmannAcc,[t2(end) t2(end)+T/5], [y2(end,1) y2(end,2) y2(end,3) y2(end,4)+xDV y2(end,5)+yDV y2(end,6)+zDV y2(end,7) y2(end,8) y2(end,9) y2(end,10) y2(end,11) y2(end,12)],options); %orbit after DV2
[X,Y,Z] = sphere(10);



% figure 
% hold on
% plot(ytotal(:,1)-ytotal(:,7)+rM,ytotal(:,2)-ytotal(:,8))
% surf(X*R_Earth, Y*R_Earth, Z*R_Earth)
% %surf(X*R_Moon+rM, Y*R_Moon, Z*R_Moon)
% xlabel('x [m]')
% ylabel('y [m]')
% title('Transfer to the Moon-centered frozen orbit')
% axis equal
% hold off


figure
hold on

plot3(y3(:,1)-y3(:,7),y3(:,2)-y3(:,8),y3(:,3)-y3(:,9),'g') %Transfer


hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Transfer to the Moon-centered frozen orbit')
axis equal


%Define Event Function, target orbit at 1000 km above Moon surface
function [value,isterminal,direction] = CrossMoonOrbit(t2,y2)
value = sqrt((y2(1)-y2(7))^2 + (y2(2)-y2(8))^2 + (y2(3)-y2(9))^2)-3366000;
isterminal = 1;
direction = 0;
end
