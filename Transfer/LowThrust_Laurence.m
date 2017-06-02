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
theta_SC = 320; %deg start angle of satellite around the Moon, measured with respect to origin of the Mooon in the positive x-axis


%orbit of the Moon
rM = 385000600; %m
vM = sqrt(Mu_Earth/rM); %m/s speed of rotation of the Moon

%orbit at the high Moon
h3 = 9000000; %m
r3 = R_Moon + h3; %m
v3 = sqrt(Mu_Moon/r3); %m/s


%options of the integrator
options1 = odeset('RelTol', 1e-9);
options2 = odeset('RelTol', 1e-9, 'Events', @L1fromEarth);
options3 = odeset('RelTol', 1e-9, 'Events', @TwoEvents);

T = 5000000; %s time of simulation
Torbitmoon = 2*pi*sqrt((r3^3)/Mu_Moon)

%{
for theta_SC = linspace(0,360,37) 
    [t1,y1] = ode45(@StartingOrbitAcc,[ 0 2*Torbitmoon ],[rM+r3*cos(theta_SC*(pi/180)) r3*sin(theta_SC*(pi/180)) 0  -v3*sin(theta_SC*(pi/180)) vM+(v3*cos(theta_SC*(pi/180))) 0 rM 0 0 0 vM 0],options1);
    [t,y,tend2] = ode45(@LowThrustAccRev,[0 -T*2],[y1(1,1) y1(1,2) y1(1,3) y1(1,4) y1(1,5) y1(1,6) y1(1,7) y1(1,8) y1(1,9) y1(1,10) y1(1,11) y1(1,12)],options3);
    if sqrt(y(1))^2 + (y(2)-y(8))^2 + (y(3)-y(9))^2)-1737000;)
end
 %}

[t1,y1] = ode45(@StartingOrbitAcc,[ 0 2*Torbitmoon ],[rM+r3*cos(theta_SC*(pi/180)) r3*sin(theta_SC*(pi/180)) 0  -v3*sin(theta_SC*(pi/180)) vM+(v3*cos(theta_SC*(pi/180))) 0 rM 0 0 0 vM 0],options1);
[t,y,tend2] = ode45(@LowThrustAccRev,[0 -T*2],[y1(1,1) y1(1,2) y1(1,3) y1(1,4) y1(1,5) y1(1,6) y1(1,7) y1(1,8) y1(1,9) y1(1,10) y1(1,11) y1(1,12)],options3);    
%[t,y,tend2] = ode45(@LowThrustAccRev,[0 -T*1.5],[r3*cos(theta_SC*(pi/180)) rM+r3*sin(theta_SC*(pi/180)) 0  -vM-v3*sin(theta*(pi/180)) v3*cos(theta*(pi/180)) 0 0 rM 0 -vM 0 0],options3);

%%%%%% Calculate Transfer times
%tyears= tend1/(3600*24*365); %calculate transfer time in years
%t2years = tend2/(3600*24*365);



%defining the sphere of the Earth
[X,Y,Z] = sphere(40);

%plotting the figure
figure
hold on
axis equal


%%%% plotting in Moon reference frame %%%%
%{
plot(y1(:,1)-y1(:,7), y1(:,2)-y1(:,8),'r')
plot(y(:,1)-y(:,7), y(:,2)-y(:,8),'b')
plot(y1(1,1)-y1(1,7), y1(1,2)-y1(1,8),'g*','DisplayName','Entering Circular Orbit')
%}
%plot3(y(end,1),y(end,2),y(end,3),'b*','DisplayName','start location S/C') %mark begin location of the satellite

%
%plot3(y(1,1),y(1,2),y(1,3),'bh','DisplayName','last location S/C') %mark last location of the satellite
plot3(y(end,1),y(end,2),y(end,3),'b*','DisplayName','start location S/C') %mark begin location of the satellite
plot3(y(:,7),y(:,8),y(:,9),'r','DisplayName','orbit Moon') %orbit of the Moon 
plot3(y1(end,7),y1(end,8),y1(end,9),'ro','DisplayName','end location Moon') %mark starting location of the Moon
plot3(y(end,7),y(end,8),y(end,9),'r*','DisplayName','start location Moon') %mark last location of the Moon
plot3(y(:,1),y(:,2),y(:,3),'b','DisplayName','Transfer Orbit') %Transfer orbit
plot3(y1(:,1),y1(:,2),y1(:,3),'g','DisplayName','Circular Moon Orbit') %Transfer orbit
plot3(y1(1,1),y1(1,2), y1(1,3),'go','DisplayName','Entering Circular Orbit')
%}

surf(X*R_Earth,Y*R_Earth, Z*R_Earth)
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
L1fromMoon = sqrt((y(1)-y(7))^2 + (y(2)-y(8))^2 + (y(3)-y(9))^2)-326390000;
%L1fromMoon = sqrt((y(1)-y(7))^2 + (y(2)-y(8))^2 + (y(3)-y(9))^2)-58010000;
value = min(CrashMoon,L1fromMoon);
isterminal = 1;
direction = 0;
end



