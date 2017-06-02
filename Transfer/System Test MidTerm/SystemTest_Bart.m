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

%Target Orbit
h2 = 10000000;
V2 = sqrt(Mu_Earth/(R_Earth+h2));
r2 = R_Earth + h2;
T2 = 2*pi*sqrt((r2^3)/Mu_Earth);

%transfer orbit
DV1 = (sqrt(((2*Mu_Earth)/r1)-(2*Mu_Earth/(r1+r2)))-sqrt(Mu_Earth/r1));
DV2 = -(sqrt(((2*Mu_Earth)/r2)-(2*Mu_Earth/(r1+r2)))-sqrt(Mu_Earth/r2));
 

options1 = odeset('RelTol', 1e-12);
options2 = odeset('RelTol',1e-12,'Events', @Time);
[t,y] = ode45(@SystemTestAcc,[0 T1],[0 -r1 0 V0 0 0 ],options1); %initial orbit
[t2,y2] = ode45(@SystemTestAcc,[t(end) t(end)+5*T2],[y(end,1) y(end,2) y(end,3) y(end,4)+DV1 y(end,5) y(end,6)],options2); 
[t3,y3] = ode45(@SystemTestAcc,[t2(end) t2(end)+T2],[y2(end,1) y2(end,2) y2(end,3) y2(end,4)-DV2 y2(end,5) y2(end,6)],options1);
time = [t; t2; t3];
ytotal = [y; y2; y3];
figure
% subplot(1,3,1)
% hold on
% plot3(y(:,1),y(:,2),y(:,3),'b') %initial orbit
% plot3(y2(:,1),y2(:,2),y2(:,3),'b') %initial orbit
% plot3(y3(:,1),y3(:,2),y3(:,3),'b')
% plot3(y2(1,1),y2(1,2),y2(1,3),'k>')
% plot3(y2(end,1),y2(end,2),y2(end,3),'k<')
% axis equal
% xlabel('x [m]')
% ylabel('y [m]')
% title('Hohmann Transfer Trajectory')
% hold off

subplot(1,2,1)
plot(time, sqrt(ytotal(:,1).^2 + ytotal(:,2).^2))
xlabel('Time [s]')
ylabel('Radius [m]')
title('Radius During Transfer')
axis([0 time(end) 0.6*10^7 1.8*10^7]) 

subplot(1,2,2)
plot(time, sqrt(ytotal(:,4).^2 + ytotal(:,5).^2))
xlabel('Time [s]')
ylabel('Velocity [m/s]')
title('Velocity During Transfer')
axis([0 time(end) 3000 10000])



%saveas(gcf,'SystemTestSimulation', 'epsc')

function [value,isterminal,direction] = Time(t,y)
value = y(1);
isterminal = 1;
direction = -1;
end


