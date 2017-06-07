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
xDV = -986.760617456355;
yDV = -400.632646774832;
zDV = 868.158790876047;
DV2 = sqrt(xDV^2 + yDV^2 + zDV^2);

%first RAAN change
DV3 = 472.8822;
xDV1 = 0.6174*DV3;
yDV1 = 0.3946*DV3;
zDV1 = -0.6806*DV3;

%second RAAN change
DV4 = 927.436;
xDV2 = 0.1611*DV4;
yDV2 = 0.7024*DV4;
zDV2 = -0.6933*DV4;

%third RAAN change
DV5 = 927.436;
xDV3 = -0.5277*DV5;
yDV3 = 0.4907*DV5;
zDV3 = -0.6933*DV5;

%fourth RAAN change
DV6 = 927.436;
xDV4 = 0.6889*DV6;
yDV4 = 0.2117*DV6;
zDV4 = 0.6933*DV6;

%options of the integrator
options = odeset('RelTol', 1e-18);
options1 = odeset('RelTol', 1e-18, 'Events', @CrossMoonOrbit);
options2 = odeset('RelTol', 1e-18, 'Events', @FirstRAAN);
options3 = odeset('RelTol', 1e-18, 'Events', @SecondRAAN);
options4 = odeset('RelTol', 1e-18, 'Events', @ThirdRAAN);
options5 = odeset('RelTol', 1e-18, 'Events', @FourthRAAN);
options6 = odeset('RelTol', 1e-18, 'Events', @RAAN);

[t,y] = ode113(@HohmannAcc,[0 T/5],[0 r1 0 -V0 0 0 rM 0 0 0 vM 0],options); %initial orbit
[t2,y2,te] = ode113(@HohmannAcc,[0 T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V1*sin(theta*(pi/180)) V1*cos(theta*(pi/180)) 0 rM 0 0 0 vM 0],options1); %Transfer Orbit
[t3,y3,te2] = ode113(@HohmannAcc,[t2(end) t2(end)+T], [y2(end,1) y2(end,2) y2(end,3) y2(end,4)+xDV y2(end,5)+yDV y2(end,6)+zDV y2(end,7) y2(end,8) y2(end,9) y2(end,10) y2(end,11) y2(end,12)],options2); %orbit after DV2
[t4,y4] = ode113(@HohmannAcc, [t3(end) t3(end)+T], [y3(end,1) y3(end,2) y3(end,3) y3(end,4)+xDV1 y3(end,5)+yDV1 y3(end,6)+zDV1 y3(end,7) y3(end,8) y3(end,9) y3(end,10) y3(end,11) y3(end,12)],options3);
[t5,y5] = ode113(@HohmannAcc, [t4(end) t4(end)+T], y4(end,:),options3);
[t6,y6] = ode113(@HohmannAcc, [t5(end) t5(end)+T], [y4(end,1) y4(end,2) y4(end,3) y4(end,4)+xDV2 y4(end,5)+yDV2 y4(end,6)+zDV2 y4(end,7) y4(end,8) y4(end,9) y4(end,10) y4(end,11) y4(end,12)],options4);
[t7,y7] = ode113(@HohmannAcc, [t6(end) t6(end)+T], [y6(end,1) y6(end,2) y6(end,3) y6(end,4)+xDV3 y6(end,5)+yDV3 y6(end,6)+zDV3 y6(end,7) y6(end,8) y6(end,9) y6(end,10) y6(end,11) y6(end,12)],options5);
[t8,y8] = ode113(@HohmannAcc, [t7(end) t7(end)+T/5],[y7(end,1) y7(end,2) y7(end,3) y7(end,4)+xDV4 y7(end,5)+yDV4 y7(end,6)+zDV4 y7(end,7) y7(end,8) y7(end,9) y7(end,10) y7(end,11) y7(end,12)],options);

[X,Y,Z] = sphere(40);

figure
hold on
plot3(y(:,1),y(:,2),y(:,3),'b')
plot3(y2(:,1),y2(:,2),y2(:,3),'b')
plot3(y3(:,1),y3(:,2),y3(:,3),'g')
plot3(y3(:,7),y3(:,8),y3(:,9),'r')
plot3(y2(:,7),y2(:,8),y2(:,9),'r')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Transfer to the Moon-centered frozen orbit')
axis equal

figure
hold on
plot3(y3(:,1)-y3(:,7),y3(:,2)-y3(:,8),y3(:,3)-y3(:,9))
plot3(y4(:,1)-y4(:,7),y4(:,2)-y4(:,8),y4(:,3)-y4(:,9))
plot3(y5(:,1)-y5(:,7),y5(:,2)-y5(:,8),y5(:,3)-y5(:,9))
plot3(y6(:,1)-y6(:,7),y6(:,2)-y6(:,8),y6(:,3)-y6(:,9))
plot3(y7(:,1)-y7(:,7),y7(:,2)-y7(:,8),y7(:,3)-y7(:,9))
plot3(y8(:,1)-y8(:,7),y8(:,2)-y8(:,8),y8(:,3)-y8(:,9))
plot3(-2239538.05808371,-672206.984197519,-2421273.02400838,'c*')
plot3(-521800,-2274600,-2425700,'*c')
plot3(1709000,-1589200,-2425700,'*c')
plot3(-2230700,-685400,2425700,'*c')
surf(X*R_Moon,Y*R_Moon,Z*R_Moon)
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal
axis vis3d
grid on

ytotal = [y3;y4;y5;y6;y7;y8];
ttotal = [t3;t4;t5;t6;t7;t8];

rx = ytotal(:,1)-ytotal(:,7);
ry = ytotal(:,2)-ytotal(:,8);
rz = ytotal(:,3)-ytotal(:,9);
vx = ytotal(:,4)-ytotal(:,10);
vy = ytotal(:,5)-ytotal(:,11);
vz = ytotal(:,6)-ytotal(:,12);

h1 = (ry.*vz)-(rz.*vy);
h2 = (rz.*vx)-(rx.*vz);
h3 = (rx.*vy)-(ry.*vx);
i = (acos(h3./(sqrt(h1.^2 + h2.^2 + h3.^2))))*(180/pi);

figure
plot(ttotal,i)
xlabel('time [s]')
ylabel('inclination [deg]')
title('Evolution of the inclination during the manoevres')

rabs = sqrt(rx.^2 + ry.^2 + rz.^2);
vabs = sqrt(vx.^2 + vy.^2 + vz.^2);
e1 = (((vy.*h3)-(vz.*h2))/Mu_Moon)-(rx./rabs);
e2 = (((vz.*h1)-(vx.*h3))/Mu_Moon)-(ry./rabs);
e3 = (((vx.*h2)-(vy.*h1))/Mu_Moon)-(rz./rabs);
e = sqrt(e1.^2 + e2.^2 + e3.^2);

figure
plot(ttotal,e)
xlabel('time [s]')
ylabel('eccentricity [-]')
title('Evolution of the eccentricity during the manoeuvres')

a = 1/((2*(rabs.^(-1)))-(vabs.*(1/Mu_Moon)));

figure
plot(ttotal,rabs)
xlabel('time [s]')
ylabel('semi-major axis')
title('Evolution of the semi-major axis during the manoeuvres')

n1 = -h2.*1;
n2 = h1.*1;
n3 = 0;
raan = [];
for i = linspace(1,length(n1),length(n1))
    if n2(i)<0
        RAAN1 = 2*pi - (acos(n1(i)/sqrt(n1(i)^2 + n2(i)^2)));
    else
        RAAN1 = acos(n1(i)/sqrt(n1(i)^2 + n2(i)^2));
    end
    if RAAN1>pi
        RAAN1 = RAAN1-(2*pi);
    end
    raan = [raan;RAAN1];
end

raan = raan*(180/pi);



figure
plot(ttotal,raan)
xlabel('time [s]')
ylabel('RAAN [deg]')
title('Evolution of the RAAN during the manoeuvres')

%Define Event Function, target orbit at 1000 km above Moon surface
function [value,isterminal,direction] = CrossMoonOrbit(t2,y2)
value = sqrt((y2(1)-y2(7))^2 + (y2(2)-y2(8))^2 + (y2(3)-y2(9))^2)-3366000;
isterminal = 1;
direction = 0;
end

function [value,isterminal,direction] = FirstRAAN(t,y)
value = sqrt((y(1)-y(7)+2239538.05808371)^2 +(y(2)-y(8)+672206.984197519)^2 + (y(3)-y(9)+2421273.02400838)^2)-37000;
isterminal = 1;
direction = -1;
end

function [value,isterminal,direction] = SecondRAAN(t,y)
value = sqrt((y(1)-y(7)+521800)^2 +(y(2)-y(8)+2274600)^2 + (y(3)-y(9)+2425700)^2)-100000;
isterminal = 1;
direction = -1;
end

function [value,isterminal,direction] = ThirdRAAN(t,y)
value = sqrt((y(1)-y(7)-1709000)^2 +(y(2)-y(8)+1589200)^2 + (y(3)-y(9)+2425700)^2)-63460;
isterminal = 1;
direction = -1;
end

function [value,isterminal,direction] = FourthRAAN(t,y)
value = sqrt((y(1)-y(7)+2230700)^2 +(y(2)-y(8)+685400)^2 + (y(3)-y(9)-2425700)^2)-184680;
isterminal = 1;
direction = -1;
end

function [value,isterminal,direction] = RAAN(t,y)
value = y(3)-y(9);
isterminal=1;
direction=1;
end
