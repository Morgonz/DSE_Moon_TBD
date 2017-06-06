clear all;
close all;
clc;

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Sun = 1.327124E+20; %m^3 s^-2 from wikipedia
Mu_Moon = 4.902801e12; %m^3 s^-2
R_Earth = 6371000; %m
R_Sun = 695508000; %m https://solarsystem.nasa.gov/planets/sun/facts
R_Moon = 1737000; %m



%Circular Orbit initial conditions
h = 500000; %m
vS = sqrt(Mu_Earth/(R_Earth + h)); %m/s
rS = R_Earth + h; %m

%Earth initial conditions
rE = 149.60E+9; %m radius of the Earth https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
vE = sqrt((Mu_Earth+Mu_Sun)/rE); %rotational speed of the Earth

%Moon initial conditions
rM = 385000600; %m
vM = sqrt((Mu_Earth+Mu_Moon)/rM); %m/s speed of rotation of the Moon

%options of the integrator
options = odeset('RelTol', 1e-12);
T = 10000000; %s

%integrator
[t,y] = ode113(@SunEarthAcc, [0 T],[rE 0 0 0 vE 0 rE+rM 0 0 0 vE+vM 0 rE+rS 0 0 0 vE+vS 0],options);


%Earth-Sun Lagrangian points
RSE = sqrt(y(:,1).^2 + y(:,2).^2 + y(:,3).^2);
rLS = RSE.*(Mu_Earth/(3*Mu_Sun))^(1/3);
xLS1 = y(:,1)-(rLS.*(y(:,1)./RSE));
yLS1 = y(:,2)-(rLS.*(y(:,2)./RSE));
zLS1 = y(:,3)-(rLS.*(y(:,3)./RSE));
rLS1 = [xLS1 yLS1 zLS1 xLS1 yLS1 zLS1];
xLS2 = y(:,1)+(rLS.*(y(:,1)./RSE));
yLS2 = y(:,2)+(rLS.*(y(:,2)./RSE));
zLS2 = y(:,3)+(rLS.*(y(:,3)./RSE));
rLS2 = [xLS2 yLS2 zLS2 xLS2 yLS2 zLS2];

%Earth-Moon Lagrangian points


% figure
% hold on
% 
% plot(y(:,7)-y(:,1),y(:,8)-y(:,2))
% plot(y(:,13)-y(:,1),y(:,14)-y(:,2))
% plot(rLS1(:,1)-y(:,1),rLS1(:,2)-y(:,2))
% hold off
% axis equal

yrot = RotatingFrame(y);
figure
hold on
plot(yrot(:,1),yrot(:,2))
plot(yrot(:,7),yrot(:,8))
plot(yrot(:,13),yrot(:,14))
hold off
axis equal
