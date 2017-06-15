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

% MANIFOLDS = [];
% 
% %options of the integrator
% options = odeset('RelTol', 2.22045e-14,'Events',@Cross);
% T = 10000000; %s
% for DVx = 97.7899
%     for DVy = 108.130100
%         for DVz = -1.26001
%     
%             [t,y,te] = ode113(@SunEarthAcc, [0 -5.8*T], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749+DVx vE+426+DVy 1.260005700065+DVz],options); %stable for the longest time 
% 
%             yrot = RotatingFrameSunEarth(y);
%             if yrot(end,8)<10000000
%                 hold on
%                 b1=plot3(yrot(:,7),yrot(:,8),yrot(:,9));
%                 
%                 MANIFOLDS = [MANIFOLDS; DVx DVy DVz yrot(end,8)];
%             end
%         end
%     end
% end
% b3=plot3(yrot(1,7),yrot(1,8),yrot(1,9),'or');
% b2=plot3(rE+rL,0,0,'k*');
% b4=plot3(rE,0,0,'bo');
% xlabel('x [m]')
% ylabel('y [m]')
% zlabel('z [m]')
% axis equal
% axis vis3d
% title('Halo orbit around Sun-Earth L2')
% legend([b3,b2,b4],{'Initial position','Sun-Earth L2','Earths location'}) 

%dV_kick = 1; %in delta V
T_days = 800;
T1 = 60*60*24*T_days; %s

MAN = [];

optionsL1 = odeset('RelTol', 2.22045e-14);
optionsL2 = odeset('RelTol', 2.22045e-14, 'Events', @toofar);
for dv_x =linspace(80,100,6)
    for dv_y = linspace(110,120,6)  
        for t_manouvre =linspace(1,10,10)
            t = 60*60*24*t_manouvre; %s

            [t1,y1] = ode113(@SunEarthAcc, [0 t], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],optionsL1);
            yrot1 = RotatingFrameSunEarth(y1);
           % v_abs = sqrt(yrot1(end,10)^2 + yrot1(end,11)^2);   
        %     dv_x = dV_kick*(yrot1(end,10)/v_abs);
        %     dv_y = dV_kick*(yrot1(end,11)/v_abs);

            [t2,y2] = ode113(@SunEarthAcc, [0 -100*t], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+dv_x y1(end,11)+dv_y y1(end,12)],optionsL2);
            yrot2 = RotatingFrameSunEarth(y2);

            if sqrt((y2(end,7)-y2(end,1))^2 + (y2(end,8)-y2(end,2))^2 + (y2(end,9)-y2(end,3))^2)<7000000
                hold on
%                 plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9));
                b1=plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9));
                PILS = sqrt(dv_x^2 + dv_y^2);
                MAN = [MAN; dv_x dv_y t_manouvre PILS];
            end

        end
disp('dv_x');disp(dv_x);disp('dv_y');disp(dv_y);
    end

end
%[t3,y3] = ode113(@SunEarthAcc, [0 T1/4.5], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],optionsL1);
%yrot3 = RotatingFrameSunEarth(y3);
hold on
[X,Y,Z]=sphere(40);
%surf(X*R_Earth+rE,Y*R_Earth,Z*R_Earth)
%b5=plot3(yrot3(:,7),yrot3(:,8),yrot3(:,9),'--');
%b4=plot3(yrot2(1,7),yrot2(1,8),yrot2(1,9),'ro');
b2=plot3(rE,0,0,'bo');
b3=plot3(rE+rL,0,0,'k*');
axis equal
axis vis3d
xlabel('x [m]')
ylabel('y [m]')
title('Trajectory from Earth to Sun-Earth L2 halo orbit')
legend([b1,b2,b3],{'Satellite trajectory','Earths location','Sun-Earth L2'});


r = sqrt((y2(end,1)-y2(end,7))^2 + (y2(end,2)-y2(end,8))^2 + (y2(end,3)-y2(end,9))^2);
vc = sqrt(Mu_Earth/r);
xc = y2(end,7)-y2(end,1);
yc = y2(end,8)-y2(end,2);
zc = y2(end,9)-y2(end,3);
theta = (atan(yc/xc)-pi/2);
vxc = vc*cos(theta);
vyc = vc*sin(theta);
vxreal = vxc + y2(end,4);
vyreal = vyc + y2(end,5);
DVx1 = y2(end,10)-vxreal;
DVy1 = y2(end,11)-vyreal;

DV = sqrt(DVx1^2 + DVy1^2);

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

value=[sqrt((y(7)-y(1))^2 + (y(8)-y(2))^2 + (y(9)-y(3))^2)-6371000-500000]; %xSr1-(1.495*10^11)];
%value = 1.496547398746715e+09 - abs(xSr1-1.5110E+11);


isterminal = [1];
direction = [0];
end
