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

%orbit at the Moon
h3 =6500000; %m
r3 = R_Moon + h3; %m
v3 = sqrt(Mu_Moon/r3); %m/s

%options of the integrator
options2 = odeset('RelTol', 1e-9, 'Events', @CrossMoonOrbit);
T = 10000000000; %s time of simulation
N = 0.34*10;
inwards = [];

for theta_SC = linspace(-180,180,361)
    [t,y,te] = ode45(@(t,y)LowThrustAccVarRev(t,y,N),[0 -T], [rM+(r3*cos(theta_SC*(pi/180))) r3*sin(theta_SC*(pi/180)) 0 -v3*sin(theta_SC*(pi/180)) vM+(v3*cos(theta_SC*(pi/180))) 0 rM 0 0 0 vM 0],options2);
    if sqrt(y(end,1)^2 + y(end,2)^2) < 385000600
        inwards = [inwards; theta_SC];
    end
    disp(theta_SC)
end

% theta_SC=-24;
% [t,y,te] = ode45(@(t,y)LowThrustAccVarRev(t,y,N),[0 -T], [rM+(r3*cos(theta_SC*(pi/180))) r3*sin(theta_SC*(pi/180)) 0 -v3*sin(theta_SC*(pi/180)) vM+(v3*cos(theta_SC*(pi/180))) 0 rM 0 0 0 vM 0],options2);
%         
% figure
% hold on
% plot((y(:,1)),y(:,2))
% %plot(y(1,1)-y(1,7),y(1,2)-y(1,8),'r*')
% axis equal
% xlabel('x [m]')
% ylabel('y [m]')
% title('Phase 2: Insertion into Moon orbit')
% hold off

%Define Event Function
function [value,isterminal,direction] = CrossMoonOrbit(t,y)
value = sqrt((y(1)-y(7))^2 + (y(2)-y(8))^2)-385000000+(3.2639*10^8);
isterminal = 1;
direction = 0;
end