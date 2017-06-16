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
aM = 384400000; %m
eM = 0.0554; %Eccentricity
rMa = aM*(1+eM);%Apocentre
rMp = aM*(1-eM); %Pericentre
vM = sqrt(((2*(Mu_Earth+Mu_Moon))/rMa)-((Mu_Earth+Mu_Moon)/aM)); %m/s speed of rotation of the Moon at the pericentre

%transfer orbit
aT  = (r1+rMp)/2; 
DV1 = 3059.93304058853; %First Delta V
T = 864000; %s time of simulation
V1 = V0+DV1; %velocity just after the initial burn
theta = -122.6; %deg angle of initial position of satellite with respect to positive x-axis

%insertion in the orbit
xDV = -975.303643734388;
yDV = -374.709373257525;
zDV = 868.158790876047;
DV2 = sqrt(xDV^2 + yDV^2 + zDV^2);

%first RAAN change
DV3 = 472.8822;
xDV1 = 0.6355*DV3;
yDV1 = 0.3646*DV3;
zDV1 = -0.6806*DV3;

%second RAAN change
DV4 = 927.436;
xDV2 = 0.1945*DV4;
yDV2 = 0.6939*DV4;
zDV2 = -0.6933*DV4;

%third RAAN change
DV5 = 927.436;
xDV3 = -0.5037*DV5;
yDV3 = 0.5154*DV5;
zDV3 = -0.6933*DV5;


%options of the integrator
options = odeset('RelTol', 1e-18);
options1 = odeset('RelTol', 1e-18, 'Events', @CrossMoonOrbit);
options2 = odeset('RelTol', 1e-18, 'Events', @FirstRAAN);
options3 = odeset('RelTol', 1e-18, 'Events', @SecondRAAN);
options4 = odeset('RelTol', 1e-18, 'Events', @ThirdRAAN);


[t,y] = ode113(@HohmannAcc,[0 T/5],[0 r1 0 -V0 0 0 rMa 0 0 0 vM 0],options); %initial orbit
[t2,y2,te] = ode113(@HohmannAcc,[0 100*T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V1*sin(theta*(pi/180)) V1*cos(theta*(pi/180)) 0 rMa 0 0 0 vM 0],options1); %Transfer Orbit
[t3,y3,te2] = ode113(@HohmannAcc,[t2(end) t2(end)+T], [y2(end,1) y2(end,2) y2(end,3) y2(end,4)+xDV y2(end,5)+yDV y2(end,6)+zDV y2(end,7) y2(end,8) y2(end,9) y2(end,10) y2(end,11) y2(end,12)],options2); %orbit after DV2
[t4,y4,te3] = ode113(@HohmannAcc,[t3(end) t3(end)+T], [y3(end,1) y3(end,2) y3(end,3) y3(end,4)+xDV1 y3(end,5)+yDV1 y3(end,6)+zDV1 y3(end,7) y3(end,8) y3(end,9) y3(end,10) y3(end,11) y3(end,12)],options3);
[t5,y5,te4] = ode113(@HohmannAcc,[t4(end) t4(end)+T], [y4(end,1) y4(end,2) y4(end,3) y4(end,4)+xDV2 y4(end,5)+yDV2 y4(end,6)+zDV2 y4(end,7) y4(end,8) y4(end,9) y4(end,10) y4(end,11) y4(end,12)],options4);
[t6,y6,te5] = ode113(@HohmannAcc,[t5(end) t5(end)+T/5], [y5(end,1) y5(end,2) y5(end,3) y5(end,4)+xDV3 y5(end,5)+yDV3 y5(end,6)+zDV3 y5(end,7) y5(end,8) y5(end,9) y5(end,10) y5(end,11) y5(end,12)],options);

[X,Y,Z] = sphere(40);

% yrot = RotatingFrame(y);
% y2rot = RotatingFrame(y2);
% y3rot = RotatingFrame(y3);

% figure
% hold on
% b1 = plot3(yrot(:,1),yrot(:,2),yrot(:,3),'r');
% plot3(y2rot(:,1),y2rot(:,2),y2rot(:,3),'r')
% plot3(y3rot(:,1),y3rot(:,2),y3rot(:,3),'r')
% b2 = plot3(y3rot(:,7),y3rot(:,8),y3rot(:,9),'color',0.5*[1 1 1]);
% plot3(y2rot(:,7),y2rot(:,8),y2rot(:,9),'color',0.5*[1 1 1])
% b6 = plot3(y3rot(end,7),y3rot(end,8),y3rot(end,9),'.','color',0.5*[1 1 1],'MarkerSize',10);
% b5 = plot3(y2rot(1,7),y2rot(1,8),y2rot(1,9),'o','color',0.5*[1 1 1]);
% b3 = plot3(y2rot(1,1),y2rot(1,2),y2rot(1,3),'k>');
% plot3(y2rot(end,1),y2rot(end,2),y2rot(end,3),'k<')
% b4 = plot3(0,0,0,'b.','MarkerSize',20);
% legend([b1,b2,b5,b6,b3,b4],{'Satellite trajectory','Moon trajectory','Moons starting position','Moons end position','Burns','Earth'})
% hold off
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')
% title('Transfer to the Moon-centered frozen orbit')
% axis equal

% y0 = [y3;y4;y5;y6];
% t0 = ([t3;t4;t5;t6]-t2(end))/(24*60*60);
% kepler = getkepler(y0(:,1)-y0(:,7),y0(:,2)-y0(:,8),y0(:,3)-y0(:,9),y0(:,4)-y0(:,10),y0(:,5)-y0(:,11),y0(:,6)-y0(:,12));
% 
% 
% for i = linspace(1,length(kepler.RAAN),length(kepler.RAAN))
%     if kepler.RAAN(i)>180*(pi/180)
%         kepler.RAAN(i) = kepler.RAAN(i)-360*(pi/180);
%     end
% end
% 
% 
figure
hold on
%b1=plot3(y3(:,1)-y3(:,7),y3(:,2)-y3(:,8),y3(:,3)-y3(:,9));
b2=plot3(y4(:,1)-y4(:,7),y4(:,2)-y4(:,8),y4(:,3)-y4(:,9));
b3=plot3(y5(:,1)-y5(:,7),y5(:,2)-y5(:,8),y5(:,3)-y5(:,9));
%b4=plot3(y6(:,1)-y6(:,7),y6(:,2)-y6(:,8),y6(:,3)-y6(:,9));
b5=plot3(y5(1,1)-y5(1,7),y5(1,2)-y5(1,8),y5(1,3)-y5(1,9),'k>');
b6=surf(X*R_Moon,Y*R_Moon,Z*R_Moon);
legend([b2,b3,b5,b6],{'Orbit after first plane change','Orbit after second plane change','Burn','Moon'});
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('First plane change')
axis equal
axis vis3d
grid on
% 
% 
% 
% figure
% 
% title('Change of the Kepler elements')
% subplot(2,2,1)
% plot(t0,kepler.INC*(180/pi))
% title('Change of inclination over time')
% xlabel('Time [days]')
% ylabel('Inclination [\circ]')
% hline = refline([0 46]);
% set(hline,'LineStyle',':')
% hline1 = refline([0 50.2]);
% set(hline1,'LineStyle',':')
% 
% 
% subplot(2,2,2)
% plot(t0,kepler.RAAN*(180/pi))
% title('Change of RAAN over time')
% xlabel('Time [days]')
% ylabel('RAAN [\circ]')
% hline2 = refline([0 284.3418-360]);
% set(hline2,'LineStyle',':')
% hline3 = refline([0 314.3418-360]);
% set(hline3,'LineStyle',':')
% hline4 = refline([0 14.3418]);
% set(hline4,'LineStyle',':')
% hline5 = refline([0 74.3418]);
% set(hline5,'LineStyle',':')
% 
% subplot(2,2,3)
% plot(t0,kepler.ECC)
% title('Change of eccentricity over time')
% xlabel('Time [days]')
% ylabel('Eccentricity [-]')
% 
% subplot(2,2,4)
% plot(t0,kepler.SMA)
% title('Change of semi-major axis over time')
% xlabel('Time [days]')
% ylabel('Semi-major axis [m]')
% hline6 = refline([0 3366000]);
% set(hline6,'LineStyle',':')


%Define Event Function, target orbit above Moon surface
function [value,isterminal,direction] = CrossMoonOrbit(t2,y2)
value = sqrt((y2(1)-y2(7))^2 + (y2(2)-y2(8))^2 + (y2(3)-y2(9))^2)-3366000;
isterminal = 1;
direction = 0;
end

function [value,isterminal,direction] = FirstRAAN(t,y)
value = sqrt((y(1)-y(7)+2269100)^2 +(y(2)-y(8)+564400)^2 + (y(3)-y(9)+2421300)^2)-8126.5;
isterminal = 1;
direction = -1;
end

function [value,isterminal,direction] = SecondRAAN(t,y)
value = sqrt((y(1)-y(7)+629900)^2 +(y(2)-y(8)+2247100)^2 + (y(3)-y(9)+2425700)^2)-5250.8;
isterminal = 1;
direction = -1;
end

function [value,isterminal,direction] = ThirdRAAN(t,y)
value = sqrt((y(1)-y(7)-1631100)^2 +(y(2)-y(8)+1669000)^2 + (y(3)-y(9)+2425700)^2)-11254.6;
isterminal = 1;
direction = -1;
end

