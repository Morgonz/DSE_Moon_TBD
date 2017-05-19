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
vM = sqrt(Mu_Earth/rM); %m/s speed of rotation of the Moon

%options of the integrator
options2 = odeset('RelTol', 1e-9, 'Events', @CrossMoonOrbit);
T = 10000000000; %s time of simulation
Speed = [];
N = 0.43;
% for N = linspace(0.1,1,91)
%     [t,y,te] = ode45(@(t,y)LowThrustAccVar(t,y,N),[0 T], [r1 0 0 0 V0 0 rM 0 0 0 vM 0],options2);
%     if abs(te)>10
%         Speed = [Speed; N te y(end,:)];
%     end
%     disp(N)
% end
[t,y,te] = ode45(@(t,y)LowThrustAccVar(t,y,N),[0 T], [r1 0 0 0 V0 0 rM 0 0 0 vM 0],options2);
Vstart = sqrt(y(1,4)^2 + y(1,5)^2);
Vend = sqrt(y(end,4)^2 + y(end,5)^2);
disp(Vstart-Vend)
% [X,Y,Z] = sphere(40);
% figure
% hold on
% plot(y(:,1),y(:,2))
% surf(X*R_Earth, Y*R_Earth, Z*R_Earth)
% axis equal
% xlabel('x [m]')
% ylabel('y [m]')
% title('Phase 1: Spiral out from Earth')
% hold off
      

%Define Event Function
function [value,isterminal,direction] = CrossMoonOrbit(t,y)
value = sqrt((y(1))^2 + (y(2))^2 + (y(3))^2)-(3.2639*10^8);
isterminal = 1;
direction = 0;
end