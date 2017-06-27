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

%Orbital parameters 
rHill = rM*(Mu_Moon/(3*Mu_Earth))^(1/3);
rL1 = rM-rHill;
rL2 = rM+rHill;
r_halo_orbit = 44000000
rot_speed_Moon = vM/rM;
v_L1 = rot_speed_Moon*(rL1);

%options of the integrator
options1 = odeset('RelTol', 2.22045e-14);
options2 = odeset('RelTol', 2.22045e-14, 'Events', @TwoEvents);
options3 = odeset('RelTol', 2.22045e-14, 'Events', @toofar_in_y);

%define delta v of kick adn total time
%dV_kick = 20 %in delta V
%T_days =60;
T_days = 65;
T = 60*60*24*T_days; %s

%integrator   (for t_manouvre = linspace(3.7,3.9,21)) 
%% The big loop
%for dV_kick = linspace(10,100,10)
dV_kick = 10
    %for t_manouvre = linspace(0.1,12.1,25)
    %for t_manouvre = linspace(0.000001,12.5,26)
    %for t_manouvre = [0.000001 0.5 1 2  4  6  8  10 11 12] 
    %for t_manouvre = [0.51  0.525   0.540] 
    for t_manouvre = [0.000001 0.525 10.3 11 12]
        %t_manouvre = 8.5
        t = 60*60*24*t_manouvre; %s

        %[t1,y1] = ode113(@EarthMoonAcc,[0 -t],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1);
        %yrot1 = RotatingFrameSunEarth(y1);
        
        %mirrored orbit
        [t1,y1] = ode113(@EarthMoonAcc,[0 -t],[rM 0 0 0 vM 0  rL1 0 -r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1); %initial orbit
        yrot1 = RotatingFrameSunEarth(y1);
        
        
        % calculating velocities
        V_abs = sqrt(y1(end,10)^2+y1(end,11)^2+y1(end,11)^2);
        V_abs_xy = sqrt(y1(end,10)^2+y1(end,11)^2);
        dv_x = (y1(end,10)/V_abs)*dV_kick;
        dv_y = (y1(end,11)/V_abs)*dV_kick;
        dv_z = (y1(end,12)/V_abs)*dV_kick;

        [t2,y2] = ode113(@EarthMoonAcc, [-t -T], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+dv_x y1(end,11)+dv_y y1(end,12)+dv_z],options2);
        %only a kick in x direction
        %[t2,y2] = ode113(@EarthMoonAcc, [-t -T], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+50 y1(end,11) y1(end,12)],options2);
        yrot2 = RotatingFrameSunEarth(y2);
        
        
        only_above_radius = 4E+8;    %define start radius of possible trajectories
        
         if yrot2(end,7)> only_above_radius
        
            % Convert from inertial to polar positions for a point
            [theta_rad, rho] = cart2pol(y2(:,7),y2(:,8));
            y2_polar = y2;
            y2_polar(:,1) = theta_rad;
            y2_polar(:,2) = rho;
            start_location = min(find(rho>only_above_radius));
            y2_polar_plot = y2_polar(start_location:end,:);
       
            for i = linspace(1,length(y2_polar_plot(:,1)),length(y2_polar_plot(:,1)));
                theta =  y2_polar_plot(i,1);
                %transfrom to inertial to polar velocities
                T_rotate = [ cos(theta) -sin(theta); sin(theta) cos(theta) ];

                V_cartesian = [y2_polar_plot(i,10) ; y2_polar_plot(i,11)];
                V_polar = inv(T_rotate)*V_cartesian;
                y2_polar_plot(i,10) = -V_polar(2);
                y2_polar_plot(i,11) = V_polar(1);             
            end
            
            
            %{
            n= 50 %taking every n'th element;
            for ii = 1:n
                y2_polar_plot_row_select =  y2_polar_plot(1:n:end,:);
            end
            %}
            y2_polar_plot_row_select =  y2_polar_plot;
            
            display(['At t_manouvre of: ', + num2str(t_manouvre)])
            hold on
            figure(1)   %intertial frame
            %plot3(y1(:,7),y1(:,8),y1(:,9),'DisplayName','sat part 1');
            plot3(y2(:,7),y2(:,8),y2(:,9),'DisplayName',['t-manouvre at '  num2str(t_manouvre)]);
            hold on
            
            %%% plotting Poincaré 
            figure(2);   %dx/dt-y plot
            plot(y2_polar_plot_row_select(:,2),y2_polar_plot_row_select(:,10),'DisplayName',['t-manouvre at '  num2str(t_manouvre)]);
            %plot(y2_polar_plot_row_select(end,2),y2_polar_plot_row_select(end,10),'ro','DisplayName','part 3'); %only last point
            hold on
            figure(3);   %dy/dt-y plot
            plot(y2_polar_plot_row_select(:,2),y2_polar_plot_row_select(:,11),'DisplayName',['t-manouvre at '  num2str(t_manouvre)]);
            %plot(y2_polar_plot_row_select(end,2),y2_polar_plot_row_select(end,11),'ro','DisplayName','part3');  %only last point
            hold on

            %}
            %{
            %%% plotting Poincaré OLD
            figure(2);   %dx/dt-y plot
            plot(y2(end,8),y2(end,10),'k*');
            hold on
            figure(3);   %dy/dt-y plot
            plot(y2(end,8),y2(end,11),'ko');
            hold on
            %}
            
            %plotting rotating frame
            figure(4)   %rotating frame
            %plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9));
            plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9),'DisplayName',['t-manouvre at '  num2str(t_manouvre)]');
            hold on
         end

    end
%end
hold on

%% Plotting of inertial frame
figure(1)
[X,Y,Z] = sphere(20);
%L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
axis equal
legend('show')
%legend(Earth)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title(['Earth fixed intertial frame with manifold trajectories iniated by kick of ' num2str(dV_kick) ' deltaV'])

hold on
%% Poin caré plotting
figure(2);

%axis equal
legend('show')
xlabel('y [m]')
ylabel('dx/dt [m/s]')
title(['Poincaré plots of dx/dt-y '])


figure(3);
%axis equal
legend('show')
xlabel('y [m]')
ylabel('dy/dt [m/s]')
title(['Poincaré plots of dy/dt-y '])

hold on
%% Plotting of rotating frame
%
figure(4)
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
%}


hold off

function [value,isterminal,direction] = TwoEvents(t,y)
crashmoon = 1737400 - sqrt((y(7)-y(1))^2 + (y(8)-y(2))^2 + (y(9)-y(3))^2);
y_axis_cross = 5E+8 - sqrt((y(7)^2)+(y(8)^2));

value = [crashmoon; y_axis_cross];
isterminal = [1;1];
direction = [0;0];

end

function [value,isterminal,direction] = toofar_in_y(t,y)
y_boundary = y(8)-(4E+8)
value = y_boundary;
isterminal = 1;
direction = 0;

end


