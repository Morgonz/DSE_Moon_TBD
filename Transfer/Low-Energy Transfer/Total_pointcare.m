clear all;
close all;
clc;

% Celest_p2ial body paramet_p2ers
Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Sun = 1.327124E+20; %m^3 s^-2 from wikipedia
Mu_Moon = 4.902801e12; %m^3 s^-2
R_Eart = 6371000; %m
R_Sun = 695508000; %m ht_p2t_p2ps://solarsyst_p2em.nasa.gov/planet_p2s/sun/fact_p2s
R_Moon = 1737400; %m mean Moon radius

%{
%% Part 2 L2 to y_axis
%Eart_p2h init_p2ial condit_p2ions
rE = 149.60E+9; %m radius of t_p2he Eart_p2h ht_p2t_p2ps://nssdc.gsfc.nasa.gov/planet_p2ary/fact_p2sheet_p2/eart_p2hfact_p2.ht_p2ml
vE = sqrt((Mu_Earth+Mu_Sun)/rE); %rot_p2at_p2ional speed of t_p2he Eart_p2h

%Sun-Eart_p2h L2
rL = rE*(Mu_Earth/(3*Mu_Sun))^(1/3);

%opt_p2ions of t_p2he int_p2egrat_p2or
options1_p2 = odeset('Reltol', 1e-18);
options2_p2 = odeset('Reltol', 1e-18, 'Event', @toofar_p2);

%define delt_p2a v of kick adn t_p2ot_p2al t_p2ime
dV_kick_p2 = -10 %in delt_p2a V
T_days_p2 = 800;
T_p2 = 60*60*24*T_days_p2; %s

%% Loop of Part 2
for t_p2_manouvre = linspace(1,180,180)
    t_p2 = 60*60*24*t_p2_manouvre; %s

    [t1_p2,y1_p2] = ode113(@SunEarthAcc, [0 t_p2], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options1_p2);
    yrot1_p2 = RotatingFrameSunEarth(y1_p2);
    % obt_p2ain end velocit_p2y
    V_abs_p2 = sqrt(y1_p2(end,10)^2+y1_p2(end,11)^2)
    dv_x_p2 = (y1_p2(end,10)/V_abs_p2)*dV_kick_p2
    dv_y_p2 = (y1_p2(end,11)/V_abs_p2)*dV_kick_p2
    
    [t2_p2,y2_p2] = ode113(@SunEarthAcc, [0 T_p2-t_p2], [y1_p2(end,1) y1_p2(end,2) y1_p2(end,3) y1_p2(end,4) y1_p2(end,5) y1_p2(end,6) y1_p2(end,7) y1_p2(end,8) y1_p2(end,9) y1_p2(end,10)+dv_x_p2 y1_p2(end,11)+dv_y_p2 y1_p2(end,12)],options2_p2);
    yrot2_p2 = RotatingFrameSunEarth(y2_p2);
    
    if yrot2_p2(end,7)<1.514E+11 & yrot2_p2(end,8)>1E+8
        display(t_p2_manouvre)
        hold on
        %%% plotting in rotating frame
        %plot3(yrot1_p2(:,7),yrot1_p2(:,8),yrot1_p2(:,9));
        %plot3(yrot2_p2(:,7),yrot2_p2(:,8),yrot2_p2(:,9));
         
        %%% plotting Poincaré 
        figure(1);   %dx/dt-y plot
        plot(yrot2_p2(end,8),yrot2_p2(end,10),'k*');
        hold on
        figure(2);   %dy/dt-y plot
        plot(yrot2_p2(end,8),yrot2_p2(end,11),'ko');
        hold on
    end
end
%}


%% Part 3 y_axis to L1 orbit


%Moon initial conditions
rM = 385000600; %m
vM = sqrt((Mu_Earth+Mu_Moon)/rM); %m/s speed of rotation of the Moon

%Orbital parameters 
rHill = rM*(Mu_Moon/(3*Mu_Earth))^(1/3);
rL1 = rM-rHill
rL2 = rM+rHill
r_halo_orbit = 38000000
rot_speed_Moon = vM/rM
v_L1 = rot_speed_Moon*(rM-rL1)

%options of the integrator
options1_p3 = odeset('RelTol', 1e-18);
options2_p3 = odeset('RelTol', 1e-18, 'Events', @TwoEvents_p3);
options3_p3 = odeset('RelTol', 1e-18, 'Events', @toofar_in_y_p3);

%define delta v of kick adn total time
dV_kick_p3 = 10 %in delta V
T_days_p3 = 60;
T_p3 = 60*60*24*T_days_p3; %s

%% The big loop
for t_manouvre = linspace(0.1,12.1,121)
    t = 60*60*24*t_manouvre; %s
    
    [t1,y1] = ode113(@EarthMoonAcc,[0 -t],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+896.5178218 0 ],options1); %initial orbit
    yrot1 = RotatingFrameSunEarth(y1);
    % calculating velocities
    V_abs = sqrt(y1(end,10)^2+y1(end,11)^2+y1(end,11)^2)
    V_abs_xy = sqrt(y1(end,10)^2+y1(end,11)^2)
    dv_x = (y1(end,10)/V_abs)*dV_kick
    dv_y = (y1(end,11)/V_abs)*dV_kick
    dv_z = (y1(end,12)/V_abs)*dV_kick
    
    [t2,y2] = ode113(@EarthMoonAcc, [-t -T], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+dv_x y1(end,11)+dv_y y1(end,12)+dv_z],options2);
    yrot2 = RotatingFrameSunEarth(y2);
    
    if yrot2(end,7)>rM*1.1
        display(['At t_manouvre of: ', + num2str(t_manouvre)])
        hold on
        figure(1)   %intertial frame
        %plot3(y1(:,7),y1(:,8),y1(:,9),'DisplayName','sat part 1');
        plot3(y2(:,7),y2(:,8),y2(:,9),'DisplayName','sattelite trajectory');
        hold on
        
        %%% plotting Poincaré 
        figure(2);   %dx/dt-y plot
        plot(y2(end,8),y2(end,10),'r*''DisplayName','part 3');
        hold on
        figure(3);   %dy/dt-y plot
        plot(y2(end,8),y2(end,11),'ro''DisplayName','part3');
        hold on
        
        %figure(2)   %rotating frame
        %plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9));
        %plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9),'DisplayName','sattelite trajectory');
        
        hold on
        
    end

end
hold on



%% Poin caré plotting
%
figure(1)
%axis equal
%legend('show')
xlabel('y [m]')
ylabel('dx/dt [m/s]')
title(['Poincaré plots of dx/dt-y '])

figure(2)
%axis equal
%legend('show')
xlabel('y [m]')
ylabel('dy/dt [m/s]')
title(['Poincaré plots of dx/dt-y '])
%}








function [value,isterminal,direction] = toofar_p2(t,y)
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
value = 1.496547398746715e+09 - abs(xSr1-1.5110E+11);   %rL - abs(xSr1 - (rE+rL)


isterminal = 1;
direction = 0;
end




