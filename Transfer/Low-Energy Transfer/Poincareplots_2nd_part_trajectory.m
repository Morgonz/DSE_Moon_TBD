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

%Sun-Earth L2
rL = rE*(Mu_Earth/(3*Mu_Sun))^(1/3);

% SUN
[X,Y,Z] = sphere(40);

%options of the integrator
options1 = odeset('RelTol', 1e-18);
options2 = odeset('RelTol', 1e-18, 'Events', @toofar);

%define delta v of kick adn total time
dV_kick = -10 %in delta V
T_days = 800;
T = 60*60*24*T_days; %s

%%
%for dV_kick = linspace(-100,-10,10)
    for t_manouvre = linspace(1,180,1)
        t_manouvre = 150
        t = 60*60*24*t_manouvre; %s

        [t1,y1] = ode113(@SunEarthAcc, [0 t], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options1);
        yrot1 = RotatingFrameSunEarth(y1);
        % obtain end velocity
        V_abs = sqrt(y1(end,10)^2+y1(end,11)^2)
        dv_x = (y1(end,10)/V_abs)*dV_kick
        dv_y = (y1(end,11)/V_abs)*dV_kick

        [t2,y2] = ode113(@SunEarthAcc, [0 T-t], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+dv_x y1(end,11)+dv_y y1(end,12)],options2);
        yrot2 = RotatingFrameSunEarth(y2);

        if yrot2(end,7)<1.514E+11 & yrot2(end,8)>1E+8
            display(t_manouvre)
            hold on
            %%% plotting in rotation frame
            figure(1);
            plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9));
            plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9));
            hold on
            %%% plotting Poincaré 
            figure(2);   %dx/dt-y plot
            plot(yrot2(end,8),yrot2(end,10),'k*');
            hold on
            figure(3);   %dy/dt-y plot
            plot(yrot2(end,8),yrot2(end,11),'ko');
            hold on
            
        end
    end
%end
hold off

%% Rotating Frame plotting
%{
figure(1);
title(['Asymptotic trajectories leaving Earth-Sun L2 with an deltaV kick of: ' num2str(dV_kick) 'm/s'])

L2 = plot3(rE+rL,0,0,'k*','DisplayName','Earth-Sun L2')
Earth = plot3(rE,0,0,'bo','DisplayName','Earth')
%surf(X*R_Sun,Y*R_Sun, Z*R_Sun,'DisplayName','Sun position')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend([L2 Earth],{'Earth-Sun L2','Earth'})
axis equal
axis vis3d


%}
%% Poin caré plotting
figure(2);

%axis equal
%legend('show')
xlabel('y [m]')
ylabel('dx/dt [m/s]')
title(['Poincaré plots of dx/dt-y '])


figure(3);
%axis equal
%legend('show')
xlabel('y [m]')
ylabel('dy/dt [m/s]')
title(['Poincaré plots of dy/dt-y '])


%% Event function
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
%value = 1.496547398746715e+09 - abs(xSr1-1.5110E+11);   %rL - abs(xSr1 - (rE+rL)
value = xSr1 - 1.4960E+11

isterminal = 1;
direction = 0;
end


