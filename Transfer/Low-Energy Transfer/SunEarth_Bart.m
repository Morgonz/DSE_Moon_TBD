clear all;
close all;
clc;

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Sun = 1.327124E+20; %m^3 s^-2 from wikipedia
R_Earth = 6371000; %m
R_Sun = 695508000; %m https://solarsystem.nasa.gov/planets/sun/facts

%Earth initial conditions
rE = 149.60E+9; %m radius of the Earth https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
vE = sqrt((Mu_Earth+Mu_Sun)/rE); %rotational speed of the Earth

%Satellite initial conditions
h = 500000; %m
vS = sqrt(Mu_Earth/(R_Earth + h)); %m/s
rS = R_Earth + h; %m

%Sun-Earth L2
rL = rE*(Mu_Earth/(3*Mu_Sun))^(1/3);

%options of the integrator
options = odeset('RelTol', 1e-18);
T = 10000000; %s

%[t,y] = ode113(@SunEarthAcc, [0 4*T], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options); %stable for the longest time 
[t,y] = ode113(@SunEarthAcc, [0 -5*T], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 0],options);
%[t,y] = ode113(@SunEarthAcc, [0 7*T], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.0663 vE+426 1.260005700065],options); %begin and end very close together after one revolution

yrot = RotatingFrameSunEarth(y);

figure
hold on
plot3(yrot(:,7),yrot(:,8),yrot(:,9))
plot3(rE+rL,0,0,'c*')
plot3(rE,0,0,'bo')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal
axis vis3d

function [value,isterminal,direction] = Cross(t,y)
if y(1)>0
    tet = atan(y(2)/y(1));
elseif y(1)<0
    tet = atan(y(2)/y(1))+pi;
else
    if y(2)>0
        tet = pi/2;
    else
        tet = -pi/2;   
    end
end
xEr1 = y(1)*cos(-tet)-y(2)*sin(-tet); %x-coordinate of the Earth in rotating frame
yEr1 = y(1)*sin(-tet)+y(2)*cos(-tet); %y-coordinate of the Earth in rotating frame
xSr1 = y(7)*cos(-tet)-y(8)*sin(-tet); %x-coordinate of the satellite in rotating frame
ySr1 = y(7)*sin(-tet)+y(8)*cos(-tet); %y-coordinate of the satellite in rotating frame

value = ySr1;
isterminal=1;
direction=1;
end


