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

MANIFOLDS = [];

%options of the integrator
options = odeset('RelTol', 2.22045e-14,'Events',@Cross);
T = 10000000; %s
for DVx = linspace(0,100,6)
    for DVy = linspace(0,100,6)
        for DVz = -1.260005700065
    
            [t,y] = ode113(@SunEarthAcc, [0 -4*T], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749+DVx vE+426+DVy 1.260005700065+DVz],options); %stable for the longest time 

            yrot = RotatingFrameSunEarth(y);
            if yrot(end,8)<10000000
                hold on
                plot3(yrot(:,7),yrot(:,8),yrot(:,9))
                MANIFOLDS = [MANIFOLDS; DVx DVy DVz yrot(end,8)];
            end
        end
    end
end
plot3(rE+rL,0,0,'c*')
plot3(rE,0,0,'bo')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal
axis vis3d


dV_kick = 1; %in delta V
T_days = 800;
T1 = 60*60*24*T_days; %s

optionsL1 = odeset('RelTol', 2.22045e-14);
optionsL2 = odeset('RelTol', 2.22045e-14, 'Events', @toofar);

for t_manouvre = linspace(1,180,36)
    t = 60*60*24*t_manouvre; %s

    [t1,y1] = ode113(@SunEarthAcc, [0 t], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],optionsL1);
    yrot1 = RotatingFrameSunEarth(y1);
    angle = atan2(yrot1(end,8),(yrot1(end,7)-(rE+rL)));
    dv_x = dV_kick*sin(angle);
    dv_y = dV_kick*-cos(angle);

    [t2,y2] = ode113(@SunEarthAcc, [0 T1-t], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7)+dV_kick*sin(angle) y1(end,8)-dV_kick*cos(angle) y1(end,9) y1(end,10) y1(end,11) y1(end,12)],optionsL2);
    yrot2 = RotatingFrameSunEarth(y2);
    
    if yrot2(end,7)<1.514E+11
        hold on
        plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9));
        plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9));
    end
end

% r = sqrt((y(end,1)-y(end,7))^2 + (y(end,2)-y(end,8))^2 + (y(end,3)-y(end,9))^2);
% vc = sqrt(Mu_Earth/r);
% xc = y(end,7)-y(end,1);
% yc = y(end,8)-y(end,2);
% zc = y(end,9)-y(end,3);
% theta = (atan(yc/xc)-pi/2);
% vxc = -vc*cos(theta);
% vyc = -vc*sin(theta);
% vxreal = vxc + y(end,4);
% vyreal = vyc + y(end,5);
% DVx1 = y(end,10)-vxreal;
% DVy1 = y(end,11)-vyreal;
% 
% DV = sqrt(DVx1^2 + DVy1^2);
% 
% options1 = odeset('RelTol', 2.22045e-14);
% options2 = odeset('RelTol', 2.22045e-14, 'Events', @Halo);
% [tS,yS] = ode113(@SunEarthAcc, [0 T/1000],[y(end,1) y(end,2) y(end,3) y(end,4) y(end,5) y(end,6) y(end,1)+xc y(end,2)+yc y(end,3)+zc vxreal vyreal 0],options1);
% [ts1,ys1] = ode113(@SunEarthAcc, [0 T],[y(end,1) y(end,2) y(end,3) y(end,4) y(end,5) y(end,6) y(end,1)+xc y(end,2)+yc y(end,3)+zc vxreal+DVx1 vyreal+DVy1 0],options2);
% [ts2,ys2] = ode113(@SunEarthAcc, [ts1(end) ts1(end)+4*T],[ys1(end,1) ys1(end,2) ys1(end,3) ys1(end,4) ys1(end,5) ys1(end,6) ys1(end,7) ys1(end,8) ys1(end,9) ys1(end,10)-DVx ys1(end,11)-DVy ys1(end,12)-DVz],options1);
% [thalo,yhalo] = ode113(@SunEarthAcc, [0 4*T], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options); %stable for the longest time
% ySrot = RotatingFrameSunEarth(yS);
% ys1rot = RotatingFrameSunEarth(ys1);
% ys2rot = RotatingFrameSunEarth(ys2);
% yhalorot = RotatingFrameSunEarth(yhalo);
% 
% 
% plot3(ySrot(:,7),ySrot(:,8),ySrot(:,9),'r')
% plot3(ys1rot(:,7),ys1rot(:,8),ys1rot(:,9))
% %plot3(ys2rot(:,7),ys2rot(:,8),ys2rot(:,9))
% %plot3(yhalorot(:,7),yhalorot(:,8),yhalorot(:,9))
% axis equal


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

value = xSr1-xEr1;
isterminal=1;
direction=0;
end

function [value,isterminal,direction] = Close(t,y)
value = sqrt((y(7)-y(1))^2 + (y(8)-y(2))^2 + (y(9)-y(3))^2) - 7656500;
isterminal = 1;
direction = 1;
end

function [value,isterminal,direction] = Halo(t,y)
value = y(8);
isterminal = 1;
direction = 0;
end

function [value,isterminal,direction] = toofar(t,y)
%disp('value found is')
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
xEr1 = y(1)*cos(-tet)-y(2)*sin(-tet);
yEr1 = y(1)*sin(-tet)+y(2)*cos(-tet);
xSr1 = y(7)*cos(-tet)-y(8)*sin(-tet);
ySr1 = y(7)*sin(-tet)+y(8)*cos(-tet);

%value = 2*1.496547398746715e+09 - sqrt((xSr1-1.5110E+11)^2 + (ySr1)^2);
value = 1.496547398746715e+09 - abs(xSr1-1.5110E+11);


isterminal = 1;
direction = 0;
end
