clear all;
close all;
clc;

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Moon = 4.902801e12; %m^3 s^-2
R_Earth = 6371000; %m mean radius of the Earth
R_Moon = 1737400; %m mean Moon radius

%Moon initial conditions
rM = 385000600; %m
vM = sqrt((Mu_Earth+Mu_Moon)/rM); %m/s speed of rotation of the Moon

T_period_Moon_sec = 2*pi*sqrt((rM^3)/Mu_Earth);
T_period_Moon_days = T_period_Moon_sec/(3600*24)

%Earth-Moon L1, seen from the Moon
rHill = rM*(Mu_Moon/(3*Mu_Earth))^(1/3);

%Satellite intital conditions
rL1 = rM-rHill;
rL2 = rM+rHill;
r_halo_orbit = 44000000;
r_halo_orbit_old = 38000000;
rot_speed_Moon = vM/rM;
v_rM_L1 = rot_speed_Moon*(rM-rL1);
v_L1 = rot_speed_Moon*(rL1);

%options of the integrator
options1 = odeset('RelTol', 2.22045e-14);
options2 = odeset('RelTol', 2.22045e-14, 'Events', @TwoEvents);
options3 = odeset('RelTol', 2.22045e-14, 'Events', @toofar_in_y);

%define delta v of kick adn total time
dV_kick = 10 %in delta V
T_days = 50;
T = 60*60*24*T_days; %s


%[t1,y1] = ode113(@EarthMoonAcc,[0 T],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1);
%yrot1 = RotatingFrameSunEarth(y1);
%hold on
%figure(2)
%plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9),'DisplayName','halo orbit');



%integrator   (for t_manouvre = linspace(3.7,3.9,21)) 
%% The big loop
%
for t_manouvre = linspace(0.1,12.25,1)
    t_manouvre = 12
    t = 60*60*24*t_manouvre; %s
        
    %[t1,y1] = ode113(@EarthMoonAcc,[0 -t],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit_old 0 v_rM_L1+896.5178218 0 ],options1); %initial orbit
    [t1,y1] = ode113(@EarthMoonAcc,[0 -t],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1);
    yrot1 = RotatingFrameSunEarth(y1);
    hold on
    % calculating velocities
    V_abs = sqrt(y1(end,10)^2+y1(end,11)^2+y1(end,11)^2);
    V_abs_xy = sqrt(y1(end,10)^2+y1(end,11)^2);
    dv_x = (y1(end,10)/V_abs)*dV_kick;
    dv_y = (y1(end,11)/V_abs)*dV_kick;
    dv_z = (y1(end,12)/V_abs)*dV_kick;
    
    [t2,y2] = ode113(@EarthMoonAcc, [-t -T], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+dv_x y1(end,11)+dv_y y1(end,12)+dv_z],options2);
    yrot2 = RotatingFrameSunEarth(y2);
    
    %plot the intitial halo
    hold on
    figure(2)   %intertial frame
    plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9),'k','DisplayName','halo orbit');
    
    
    if yrot2(end,7)>rM*1.0
        display(['At t_manouvre of: ', + num2str(t_manouvre)])
        hold on
        figure(1)   %intertial frame
        
        plot3(y2(:,7),y2(:,8),y2(:,9),'DisplayName','sattelite trajectory');
        hold on  
        
        % Convert from inertial to polar velocities
        [theta_rad, rho] = cart2pol(y2(end,7),y2(end,8))
        T = [  cos(theta_rad) -sin(theta_rad);
               sin(theta_rad) cos(theta_rad)  ]
        V_cartesian = [y2(end,10) ; y2(end,11)]
        V_polar = inv(T)*V_cartesian
        
        
        
        figure(2)   %rotating frame
        plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9),'k');
        plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9),'DisplayName','sattelite trajectory');   
        hold on  
    end
    
    %}
end
hold on
figure(1)
Moon = plot3(y1(:,1),y1(:,2),y1(:,3), 'k' ,'DisplayName','Moon orbit');
Moon = plot3(y2(:,1),y2(:,2),y2(:,3), 'k' ,'DisplayName','Moon orbit');
Moon_end = plot3(y1(1,1),y1(1,2),y1(1,3), 'ko','DisplayName','Moon end location');
Moon_begin = plot3(y2(end,1),y2(end,2),y2(end,3), 'k*' ,'DisplayName','Moon start location');
hold on

%% Plotting of inertial frame
figure(1)
[X,Y,Z] = sphere(20);
%L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
axis equal
legend('show')
legend(Moon)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title(['Earth fixed intertial frame with manifold trajectories iniated by kick of ' num2str(dV_kick) ' deltaV'])

hold off

%}
%% Plotting of rotating frame
figure(2)
hold on
[X,Y,Z] = sphere(20);
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'k+','DisplayName','Earth-Moon L2');
%Moon = plot3(rM,0,0,'go','DisplayName','Moon');
%plot3(y(:,1),y(:,2),y(:,3));  %plotting the Moon
Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');

%legend([L1 Earth],{'Earth-Moon L1','Earth'})
axis equal
axis vis3d
legend('show')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title(['Earth fixed rotating frame with manifold trajectories initiated by kick of ' num2str(dV_kick) ' deltaV'])

hold off





function [value,isterminal,direction] = TwoEvents(t,y)
crashmoon = 1737400 - sqrt((y(7)-y(1))^2 + (y(8)-y(2))^2 + (y(9)-y(3))^2);
%{
%x_coord = y(7);
if y(8)>3.5E+8
    y_axis_cross = y(7);
else
     y_axis_cross = -10000;
end
%}
 y_axis_cross = 4E+8 - sqrt((y(7)^2)+(y(8)^2));

value = [crashmoon; y_axis_cross];
isterminal = [1;1];
direction = [0;0];

end

function [value,isterminal,direction] = toofar_in_y(t,y)
y_boundary = y(8)-(5E+8);
value = y_boundary;
isterminal = 1;
direction = 0;

end


