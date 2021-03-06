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
theta = 0; %deg angle the satellite is at the start of the simulation with respect to the positive x-axis
theta_SC = 60; %deg start angle of satellite around the Moon, measured with respect to origin of the Mooon in the positive x-axis


%orbit of the Moon
rM = 385000600; %m
vM = sqrt((Mu_Earth+Mu_Moon)/rM); %m/s speed of rotation of the Moon

%orbit at the Moon
h3 = 9000000; %m
r3 = R_Moon + h3; %m
v3 = sqrt(Mu_Moon/r3); %m/s


%options of the integrator
options1 = odeset('RelTol', 1e-9);
options2 = odeset('RelTol', 1e-9, 'Events', @L1fromEarth);
options3 = odeset('RelTol', 1e-9, 'Events', @TwoEvents);



T = 5000000; %s time of simulation

%[t,y] = ode45(@LowThrustAcc,[0 10*T],[0 r1 0 V0 0 0 0 rM 0 vM 0 0],options1); %Starting orbit of the satellite
%[t,y,tend1] = ode45(@LowThrustAcc,[0 4*T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V0*sin(theta*(pi/180)) V0*cos(theta*(pi/180))  0 rM 0 0 0 vM 0],options2); %Transfer orbit
%[t,y,tend2] = ode45(@LowThrustAccRev,[4*T 0],[r1*cos(theta*(pi/180))

[t,y,tend2] = ode45(@LowThrustAccRev,[0 -T*1.5],[r3*cos(theta_SC*(pi/180)) rM+r3*sin(theta_SC*(pi/180)) 0  -vM-v3*sin(theta*(pi/180)) v3*cos(theta*(pi/180)) 0 0 rM 0 -vM 0 0],options3);


%%%%%% Calculate Transfer times
%tyears= tend1/(3600*24*365); %calculate transfer time in years
t2years = tend2/(3600*24*365);

%defining the sphere of the Earth
[X,Y,Z] = sphere(40);

%plotting the figure
figure
%plot3(y(:,1),y(:,2),y(:,3),'b','DisplayName','S/C') %Low-Earth orbit of the satellite
hold on

%%%% plotting in Moon reference frame %%%%
%plot(y(:,1)-y(:,7), y(:,2)-y(:,8))


%
plot3(y(1,1),y(1,2),y(1,3),'bh','DisplayName','last location S/C') %mark last location of the satellite
plot3(y(end,1),y(end,2),y(end,3),'b*','DisplayName','start location S/C') %mark begin location of the satellite
plot3(y(:,7),y(:,8),y(:,9),'r','DisplayName','orbit Moon') %orbit of the Moon 
plot3(y(1,7),y(1,8),y(1,9),'ro','DisplayName','end location Moon') %mark starting location of the Moon
plot3(y(end,7),y(end,8),y(end,9),'r*','DisplayName','start location Moon') %mark last location of the Moon
plot3(y(:,1),y(:,2),y(:,3),'b') %Transfer orbit
%}

%plot3(y3(:,1),y3(:,2),y3(:,3),'c') %orbit during coasting
%plot3(y3(:,7),y3(:,8),y3(:,9),'r')%Orbit of the Moon during coasting
%surf(X*R_Earth, Y*R_Earth, Z*R_Earth) % where (a,b,c) is center of the sphere
surf(X*R_Moon, Y*R_Moon, Z*R_Moon) % where (a,b,c) is center of the sphere
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')


legend('show')
axis equal




%Define Event Function
function [value,isterminal,direction] = L1fromEarth(t,y)
value = sqrt((y(1))^2 + (y(2))^2 + (y(3))^2)-326390000;
isterminal = 1;
direction = 0;
end

%{
function [value,isterminal,direction] = L1fromMoon(t,y)
%value = sqrt((y(1)-y(7))^2 + (y(2)-y(8))^2 + (y(3)-y(9))^2)-58010000; %L1 sphere seen from the Moon
value = sqrt((y(1))^2 + (y(2))^2 + (y(3))^2)-326390000; %L1 sphere seen from the Earth
isterminal = 1;
direction = 0;
end
%}

function [value,isterminal,direction] = TwoEvents(t,y)
CrashMoon = sqrt((y(1)-y(7))^2 + (y(2)-y(8))^2 + (y(3)-y(9))^2)-1737000;
L1fromMoon = sqrt((y(1))^2 + (y(2))^2 + (y(3))^2)-326390000;
value = min(CrashMoon,L1fromMoon);
isterminal = 1;
direction = 0;
end



