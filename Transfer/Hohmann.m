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
T  = 100*2*pi*sqrt((R_Earth+h)^3/Mu_Earth); %s
r1 = R_Earth + h; %m

%orbit of the Moon
rM = 385000600; %m
vM = sqrt(Mu_Earth/rM); %m/s speed of rotation of the Moon

%transfer orbit
aT  = (r1+rM)/2;
DV1 = sqrt(Mu_Earth/r1)*(sqrt((2*rM)/(r1+rM))-1);
DV2 = sqrt(Mu_Earth/rM)*(1-sqrt((2*r1)/(r1+rM)));

%orbit around the Moon
h1 = 200000; %m
r3 = h1 + R_Moon; %radius of the orbit around the Moon
v3 = sqrt(Mu_Moon/r3); %m/s orbiting speed in Moon reference frame


%options of the integrator
options = odeset('RelTol', 1e-12);

[t,y] = ode45(@HohmannAcc,[0 T],[0 r1 0 -V0-DV1 0 0 rM 0 0 0 vM 0],options);
[t2,y2] = ode45(@HohmannAcc,[0 T],[y(length(y),1) y(length(y),2) y(length(y),3) y(length(y),4) y(length(y),5) y(length(y),6) y(length(y),7) y(length(y),8) y(length(y),9) y(length(y),10) y(length(y),11) y(length(y),12)],options);
[t3,y3] = ode45(@HohmannAcc, [0 T],[rM r3 0 -v3-DV2 vM 0 rM 0 0 0 vM 0],options);

[X,Y,Z] = sphere(40);

figure
%plot3(y(:,1),y(:,2),y(:,3),'b')
hold on
surf(X*R_Earth, Y*R_Earth, Z*R_Earth) % where (a,b,c) is center of the sphere
plot3(y(:,7),y(:,8),y(:,9),'r')
plot3(y3(:,1), y3(:,2), y3(:,3),'g')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal



