clear all;
close all;
clc;

% Celest_p2ial body paramet_p2ers
Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Sun = 1.327124E+20; %m^3 s^-2 from wikipedia
Mu_Moon = 4.902801e12; %m^3 s^-2
R_Earth = 6371000; %m
R_Sun = 695508000; %m ht_p2t_p2ps://solarsyst_p2em.nasa.gov/planet_p2s/sun/fact_p2s
R_Moon = 1737400; %m mean Moon radius

%
%% Part 2 L2 to y_axis
%Eart_p2h init_p2ial condit_p2ions
rE = 149.60E+9; %m radius of t_p2he Eart_p2h ht_p2t_p2ps://nssdc.gsfc.nasa.gov/planet_p2ary/fact_p2sheet_p2/eart_p2hfact_p2.ht_p2ml
vE = sqrt((Mu_Earth+Mu_Sun)/rE); %rot_p2at_p2ional speed of t_p2he Eart_p2h

%Sun-Eart_p2h L2
rL = rE*(Mu_Earth/(3*Mu_Sun))^(1/3);

%options of t_p2he int_p2egrat_p2or
options1_p2 = odeset('Reltol', 2.22045e-14);
options2_p2 = odeset('Reltol', 2.22045e-14, 'Event', @toofar_p2);

%define delt_p2a v of kick adn t_p2ot_p2al t_p2ime
%T_days_p2 = 600;
T_days_p2 = 400;
T_p2 = 60*60*24*T_days_p2; %s

DATA2 =[];
DATA3 =[];
yrot2_p2_plot = [];

%% Loop of Part 2
%for dV_kick_p2 = linspace(-50,10,61)
    %display(['dV_kick_p2 is ', num2str(dV_kick_p2)])
dV_kick_p2 = -10;
    for t_p2_manouvre = linspace(1,180,180)
    %for t_p2_manouvre = linspace(59.14,59.2,6)
    %for t_p2_manouvre = [ 59.14 ] 
        t_p2 = 60*60*24*t_p2_manouvre; %s

        [t1_p2,y1_p2] = ode113(@SunEarthAcc, [0 t_p2], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options1_p2);
        yrot1_p2 = RotatingFrameSunEarth(y1_p2);
        % obt_p2ain end velocit_p2y
        V_abs_p2 = sqrt(y1_p2(end,10)^2+y1_p2(end,11)^2);
        dv_x_p2 = (y1_p2(end,10)/V_abs_p2)*dV_kick_p2;
        dv_y_p2 = (y1_p2(end,11)/V_abs_p2)*dV_kick_p2;

        %[t2_p2,y2_p2] = ode113(@SunEarthAcc, [0 T_p2-t_p2], [y1_p2(end,1) y1_p2(end,2) y1_p2(end,3) y1_p2(end,4) y1_p2(end,5) y1_p2(end,6) y1_p2(end,7) y1_p2(end,8) y1_p2(end,9) y1_p2(end,10)+dv_x_p2 y1_p2(end,11)+dv_y_p2 y1_p2(end,12)],options2_p2);
        [t2_p2,y2_p2] = ode113(@SunEarthAcc, [0 T_p2], [y1_p2(end,1) y1_p2(end,2) y1_p2(end,3) y1_p2(end,4) y1_p2(end,5) y1_p2(end,6) y1_p2(end,7) y1_p2(end,8) y1_p2(end,9) y1_p2(end,10)+dv_x_p2 y1_p2(end,11)+dv_y_p2 y1_p2(end,12)],options2_p2);
            
        yrot2_p2 = RotatingFrameSunEarth(y2_p2);
        
        if yrot2_p2(end,7)<1.514E+11 & yrot2_p2(end,8)>2E+8
            
            new_data = [yrot2_p2(end,8) yrot2_p2(end,10) yrot2_p2(end,11) ];
            yrot2_p2_plot = [yrot2_p2_plot ; new_data ];
            
            T_part2 = (t2_p2(1)-t2_p2(end))/(3600*24)
            
            hold on
            %%% plotting in rotating frame
            figure(2);
            plot3(yrot1_p2(:,7),yrot1_p2(:,8),yrot1_p2(:,9));
            plot3(yrot2_p2(:,7),yrot2_p2(:,8),yrot2_p2(:,9),'DisplayName',['manouvre at', num2str(t_p2_manouvre) ]);
            
            %{
            %%% plotting Poincaré 
            figure(1);   %dx/dt-y plot
            %plot(yrot2_p2(end,8),yrot2_p2(end,10),'k*','DisplayName','part 2');
            plot(yrot2_p2(end,8),yrot2_p2(end,10),'k*');
            hold on
            figure(2);   %dy/dt-y plot
            %plot(yrot2_p2(end,8),yrot2_p2(end,11),'ko','DisplayName','part 2');
            plot(yrot2_p2(end,8),yrot2_p2(end,11),'k*');
            hold on
            %}
            
            
            
            %DATA2 = [DATA2; dV_kick_p2 t_p2_manouvre yrot2_p2(end,8) yrot2_p2(end,10) yrot2_p2(end,11)];
        end
    end
%end    
    %}

hold on
%%% plotting Poincaré 
figure(1);   %dx/dt-y plot
subplot(2,1,1)
hold on
%plot(yrot2_p2(end,8),yrot2_p2(end,10),'k*','DisplayName','part 2');
plot(yrot2_p2_plot(:,1),yrot2_p2_plot(:,2),'k*','DisplayName','dx/dt part 2 trajectory');



figure(1);   %dy/dt-y plot
hold on
subplot(2,1,2)
%plot(yrot2_p2(end,8),yrot2_p2(end,11),'ko','DisplayName','dy/dt part 2');
plot(yrot2_p2_plot(:,1),yrot2_p2_plot(:,3),'k*','DisplayName','Trajectory part 2');

    
    
    
    
%% Part 3 y_axis to L1 orbit


%Moon initial conditions
rM = 385000600; %m
vM = sqrt((Mu_Earth+Mu_Moon)/rM); %m/s speed of rotation of the Moon

%Orbital parameters 
rHill = rM*(Mu_Moon/(3*Mu_Earth))^(1/3);
rL1 = rM-rHill;
rL2 = rM+rHill;
r_halo_orbit = 44000000;
rot_speed_Moon = vM/rM;
v_L1 = rot_speed_Moon*(rL1);

%options of the integrator
options1_p3 = odeset('RelTol', 2.22045e-14);
options2_p3 = odeset('RelTol', 2.22045e-14, 'Events', @TwoEvents_p3);
options3_p3 = odeset('RelTol', 2.22045e-14, 'Events', @toofar_in_y_p3);

%define delta v of kick adn total time
%dV_kick_p3 = 10; %in delta V
T_days_p3 = 70;
T_p3 = 60*60*24*T_days_p3; %s
 
%% The big loop of part 3
%for dV_kick_p3 = linspace(-30,50,81)
dV_kick_p3 = 10;
    display(['dV_kick_p3 is ', num2str(dV_kick_p3)]);
    
    %for t_manouvre_p3 = linspace(0.000001,12.5,26)
    %for t_manouvre_p3 = linspace(0.52,0.62,6)
    %for t_manouvre_p3 = linspace(0.10,0.2,11)
    %for t_manouvre_p3 = [0.51  0.525   0.540] %linspace(0.524,0.526,3)
    %for t_manouvre_p3 = [0.000001 0.5 1 2  4  6  8  10 11 12] 
    %for t_manouvre_p3 = [0.000001 0.5 0.525 10.3 11 ]
    %for t_manouvre_p3 = [0.000001 0.525 10.3 ]
    %for t_manouvre_p3 = 0.525
    for t_manouvre_p3 = [1.4 1.5 1.55 1.6 1.65]
        t_p3 = 60*60*24*t_manouvre_p3; %s

        %[t1_p3,y1_p3] = ode113(@EarthMoonAcc,[0 -t_p3],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+896.5178218 0 ],options1_p3); %initial orbit
        [t1_p3,y1_p3] = ode113(@EarthMoonAcc,[0 -t_p3],[rM 0 0 0 vM 0  rL1 0 -r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1_p3);
        yrot1_p3 = RotatingFrameSunEarth(y1_p3);
        
        % calculating velocities
        V_abs_p3 = sqrt(y1_p3(end,10)^2+y1_p3(end,11)^2+y1_p3(end,11)^2);
        V_abs_xy_p3 = sqrt(y1_p3(end,10)^2+y1_p3(end,11)^2);
        dv_x_p3 = (y1_p3(end,10)/V_abs_p3)*dV_kick_p3;
        dv_y_p3 = (y1_p3(end,11)/V_abs_p3)*dV_kick_p3;
        dv_z_p3 = (y1_p3(end,12)/V_abs_p3)*dV_kick_p3;

        [t2_p3,y2_p3] = ode113(@EarthMoonAcc, [-t_p3 -T_p3], [y1_p3(end,1) y1_p3(end,2) y1_p3(end,3) y1_p3(end,4) y1_p3(end,5) y1_p3(end,6) y1_p3(end,7) y1_p3(end,8) y1_p3(end,9) y1_p3(end,10)+dv_x_p3 y1_p3(end,11)+dv_y_p3 y1_p3(end,12)+dv_z_p3],options2_p3);
        yrot2_p3 = RotatingFrameSunEarth(y2_p3);
        
        Time_part3 = (t2_p3(1)-t2_p3(end))/(3600*24);
        
        %only_above_radius = 4.0E+8;    %define start radius of possible trajectories
        only_above_radius = rM+10000000;
        %stop_at_radius = 3.998e+8
        %if yrot2_p3(end,7)> only_above_radius
        if sqrt((y2_p3(end,7)^2)+(y2_p3(end,8)^2)) > only_above_radius
        %if sqrt((y2_p3(end,7)^2)+(y2_p3(end,8)^2)) < stop_at_radius 
            
            % Convert from inertial to polar positions for a point
            [theta_rad, rho] = cart2pol(y2_p3(:,7),y2_p3(:,8));
            y2_polar = y2_p3;
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
            n= 5; %taking every n'th element
            for ii = 1:n
                y2_polar_plot_row_select =  y2_polar_plot(1:n:end,:);
            end
            %}
            y2_polar_plot_row_select =  y2_polar_plot;
            hold on         
            
            figure(3);   %intertial frame part 3
            hold on
            plot3(y2_p3(:,7),y2_p3(:,8),y2_p3(:,9),'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
            %Moon = plot3(y1_p3(:,1),y1_p3(:,2),y1_p3(:,3), 'k' ,'DisplayName','Moon orbit');
            %Moon = plot3(y2_p3(:,1),y2_p3(:,2),y2_p3(:,3), 'k' ,'DisplayName','Moon orbit');
            hold on 
            
            figure(4);  %rotating frame part 3
            hold on
            %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
            plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
            Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
            hold on
                        
            
            %%% plotting Poincaré 
            figure(1);   %dx/dt-y plot
            subplot(2,1,1)
            %plot(y2_p3(end,8),y2_p3(end,10),'ro','DisplayName','part 3');
            plot(y2_polar_plot_row_select(:,2),y2_polar_plot_row_select(:,10),'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
            hold on
            figure(1);   %dy/dt-y plot
            subplot(2,1,2)
            %plot(y2_p3(end,8),y2_p3(end,11),'ro','DisplayName','part3');
            plot(y2_polar_plot_row_select(:,2),y2_polar_plot_row_select(:,11),'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
            hold on
            %DATA3 = [DATA3; dV_kick_p3 t_manouvre_p3 y2_p3(end,8) y2_p3(end,10) y2_p3(end,11)];
                                  
        end

    end
%end

%hold on
%figure(4);  %rotating frame part 3
%Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
%hold on


hold off



%% Poin caré plotting
%
figure(1)
subplot(2,1,1)
%axis equal
legend('show')
xlabel('y [m]')
ylabel('dx/dt [m/s]')
axis([1E+8 5E+8 -2500 -500])
title(['Velocities in the x-direction when crossing the inertial Earth y-axis of trajectory part 2 and 3'])

figure(1)
subplot(2,1,2)
%axis equal
legend('show')
xlabel('y [m]')
ylabel('dy/dt [m/s]')
axis([1E+8 5E+8 -1500 0])
title(['Velocities in the y-direction when crossing the inertial Earth y-axis of trajectory part 2 and 3'])
%}



figure(2);
L2 = plot3(rE+rL,0,0,'k*','DisplayName','Earth-Sun L2');
Earth = plot3(rE,0,0,'bo','DisplayName','Earth');
%surf(X*R_Sun,Y*R_Sun, Z*R_Sun,'DisplayName','Sun position')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
%legend([L2 Earth],{'Earth-Sun L2','Earth'})
legend('show')
axis equal
axis vis3d

%% For plotting Moon and Sat initial
T_clean = 27.5*60*60*24;
T_clean_halo = 12.25*60*60*24;
[t_clean,y_clean] = ode113(@EarthMoonAcc,[0 -T_clean],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1_p3);
[t_clean_halo,y_clean_halo] = ode113(@EarthMoonAcc,[0 -T_clean_halo],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1_p3);

yrot_clean = RotatingFrameSunEarth(y_clean);
yrot_clean_halo = RotatingFrameSunEarth(y_clean_halo);

%%
figure(3); %inertial frame
hold on
Moon = plot3(y_clean(:,1),y_clean(:,2),y_clean(:,3), 'k--' ,'DisplayName','Moon orbit');


[X,Y,Z] = sphere(20);
%L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
axis equal
legend('show')
%legend(Moon)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

%%
figure(4); %rotating frame
[X,Y,Z] = sphere(20);
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'k+','DisplayName','Earth-Moon L2');
%Moon = plot3(rM,0,0,'go','DisplayName','Moon');
%plot3(y(:,1),y(:,2),y(:,3));  %plotting the Moon
Halo = plot3(yrot_clean_halo(:,7),yrot_clean_halo(:,8),yrot_clean_halo(:,9), 'k--' ,'DisplayName','Halo orbit around L1');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
Earth = plot3(0,0,0,'b.','MarkerSize',35,'DisplayName','Earth location');

%surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
Moon = plot3(rM,0,0,'g.','MarkerSize',25,'DisplayName','Moon location');
%legend([L1 Earth],{'Earth-Moon L1','Earth'})
axis equal
axis vis3d
legend('show')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title(['Earth fixed rotating frame with manifold trajectories initiated by kick of ' num2str(dV_kick_p3) ' deltaV'])





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


function [value,isterminal,direction] = TwoEvents_p3(t,y)
crashmoon = 1737400 - sqrt((y(7)-y(1))^2 + (y(8)-y(2))^2 + (y(9)-y(3))^2);
%x_coord = y(7);
%{
if y(8)>3.5E+8
    y_axis_cross = y(7);
else
     y_axis_cross = -10000;
end
%}
%y_axis_cross = 3.978726194287806e+08 - sqrt((y(7)^2)+(y(8)^2));

%y_axis_cross = 4.0e+08 - sqrt((y(7)^2)+(y(8)^2));
y_axis_cross = 4.5e+08 - sqrt((y(7)^2)+(y(8)^2));



value = [crashmoon; y_axis_cross];
isterminal = [1;1];
direction = [0;0];

end

function [value,isterminal,direction] = toofar_in_y_p3(t,y)
y_boundary = y(8)-(1E+8)
value = y_boundary;
isterminal = 1;
direction = 0;

end




