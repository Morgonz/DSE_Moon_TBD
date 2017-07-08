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
dV_kick_p3 = 10; %in delta V
T_days_p3 = 90;
T_p3 = 60*60*24*T_days_p3; %s

%% The big loop
%for dV_kick_p3 = linspace(-30,50,81)
dV_kick_p3 = 10;
    display(['dV_kick_p3 is ', num2str(dV_kick_p3)])
    
%% First run
%t_manouvre_p3 = [ 0.525 10.3 ]
t_manouvre_p3 = [ 0.525 ]
    t_p3 = 60*60*24*t_manouvre_p3; %s

    %[t1_p3,y1_p3] = ode113(@EarthMoonAcc,[0 -t_p3],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+896.5178218 0 ],options1_p3); %initial orbit
    [t1_p3,y1_p3] = ode113(@EarthMoonAcc,[0 -t_p3],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1_p3);
    yrot1_p3 = RotatingFrameSunEarth(y1_p3);

    % calculating velocities
    V_abs_p3 = sqrt(y1_p3(end,10)^2+y1_p3(end,11)^2+y1_p3(end,11)^2);
    V_abs_xy_p3 = sqrt(y1_p3(end,10)^2+y1_p3(end,11)^2);
    dv_x_p3 = (y1_p3(end,10)/V_abs_p3)*dV_kick_p3;
    dv_y_p3 = (y1_p3(end,11)/V_abs_p3)*dV_kick_p3;
    dv_z_p3 = (y1_p3(end,12)/V_abs_p3)*dV_kick_p3;

    [t2_p3,y2_p3] = ode113(@EarthMoonAcc, [-t_p3 -T_p3], [y1_p3(end,1) y1_p3(end,2) y1_p3(end,3) y1_p3(end,4) y1_p3(end,5) y1_p3(end,6) y1_p3(end,7) y1_p3(end,8) y1_p3(end,9) y1_p3(end,10)+dv_x_p3 y1_p3(end,11)+dv_y_p3 y1_p3(end,12)+dv_z_p3],options2_p3);
    yrot2_p3 = RotatingFrameSunEarth(y2_p3);


    figure(3);   %intertial frame part 3
    hold on
    plot3(y2_p3(:,7),y2_p3(:,8),y2_p3(:,9),'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    %Moon = plot3(y1_p3(:,1),y1_p3(:,2),y1_p3(:,3), 'k' ,'DisplayName','Moon orbit');
    %Moon = plot3(y2_p3(:,1),y2_p3(:,2),y2_p3(:,3), 'k' ,'DisplayName','Moon orbit');
    hold on 

    figure(4);  %rotating frame part 3
    subplot(2, 3, 1);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    Traj1 = plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color',[0.95 0.33 0.1],'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on
    
    figure(4);  %rotating frame part 3
    subplot(2, 3,2);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color',[0.85 0.33 0.1],'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on
    
    figure(4);  %rotating frame part 3
    subplot(2, 3,3);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color',[0.85 0.33 0.1],'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on
    
    figure(4);  %rotating frame part 3
    subplot(2, 3, 4);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color',[0.85 0.33 0.1],'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on
    
    figure(4);  %rotating frame part 3
    subplot(2, 3, 5);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color',[0.85 0.33 0.1],'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on
    
    

%% Second run

t_manouvre_p3 = [10.3 ]
    t_p3 = 60*60*24*t_manouvre_p3; %s
    
    %[t1_p3,y1_p3] = ode113(@EarthMoonAcc,[0 -t_p3],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+896.5178218 0 ],options1_p3); %initial orbit
    [t1_p3,y1_p3] = ode113(@EarthMoonAcc,[0 -t_p3],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1_p3);
    yrot1_p3 = RotatingFrameSunEarth(y1_p3);

    % calculating velocities
    V_abs_p3 = sqrt(y1_p3(end,10)^2+y1_p3(end,11)^2+y1_p3(end,11)^2);
    V_abs_xy_p3 = sqrt(y1_p3(end,10)^2+y1_p3(end,11)^2);
    dv_x_p3 = (y1_p3(end,10)/V_abs_p3)*dV_kick_p3;
    dv_y_p3 = (y1_p3(end,11)/V_abs_p3)*dV_kick_p3;
    dv_z_p3 = (y1_p3(end,12)/V_abs_p3)*dV_kick_p3;
    
    %define new time
    
    T_p3_new = (10.75)*3600*24+t_p3
    
    %[t2_p3,y2_p3] = ode113(@EarthMoonAcc, [-t_p3 -T_p3], [y1_p3(end,1) y1_p3(end,2) y1_p3(end,3) y1_p3(end,4) y1_p3(end,5) y1_p3(end,6) y1_p3(end,7) y1_p3(end,8) y1_p3(end,9) y1_p3(end,10)+dv_x_p3 y1_p3(end,11)+dv_y_p3 y1_p3(end,12)+dv_z_p3],options2_p3);
    [t2_p3_new,y2_p3] = ode113(@EarthMoonAcc, [-t_p3 -T_p3_new], [y1_p3(end,1) y1_p3(end,2) y1_p3(end,3) y1_p3(end,4) y1_p3(end,5) y1_p3(end,6) y1_p3(end,7) y1_p3(end,8) y1_p3(end,9) y1_p3(end,10)+dv_x_p3 y1_p3(end,11)+dv_y_p3 y1_p3(end,12)+dv_z_p3],options2_p3);

    yrot2_p3 = RotatingFrameSunEarth(y2_p3);

    %{
    figure(3);   %intertial frame part 3
    hold on
    plot3(y2_p3(:,7),y2_p3(:,8),y2_p3(:,9),'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    %Moon = plot3(y1_p3(:,1),y1_p3(:,2),y1_p3(:,3), 'k' ,'DisplayName','Moon orbit');
    %Moon = plot3(y2_p3(:,1),y2_p3(:,2),y2_p3(:,3), 'k' ,'DisplayName','Moon orbit');
    hold on 
    %}
    
    figure(4);  %rotating frame part 3
    subplot(2, 3, 1);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    Traj2 = plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color',[0 0.8 1] ,'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on
    
    figure(4);  %rotating frame part 3
    subplot(2, 3,2);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color', [0 0.8 1],'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on
    
    figure(4);  %rotating frame part 3
    subplot(2, 3,3);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color',[0 0.8 1],'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on
    
    figure(4);  %rotating frame part 3
    subplot(2, 3, 4);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color',[0 0.8 1],'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on
    
    figure(4);  %rotating frame part 3
    subplot(2, 3, 5);
    hold on
    %plot3(yrot1_p3(:,7),yrot1_p3(:,8),yrot1_p3(:,9),'k');
    plot3(yrot2_p3(:,7),yrot2_p3(:,8),yrot2_p3(:,9),'Color',[0 0.8 1],'DisplayName',['Trajectory part 3 with kick at '  num2str(-t_manouvre_p3) ' days'],'LineWidth',1);
    Insert_place = plot3(yrot2_p3(1,7),yrot2_p3(1,8),yrot2_p3(1,9),'ko','DisplayName','Inserting into L1 Earth-Moon orbit');   
    hold on


%% For plotting Moon and Sat initial
T_clean = 27.5*60*60*24
T_clean_halo = 12.25*60*60*24
[t_clean,y_clean] = ode113(@EarthMoonAcc,[0 -T_clean],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1_p3);
[t_clean_halo,y_clean_halo] = ode113(@EarthMoonAcc,[0 -T_clean_halo],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1_p3);

yrot_clean = RotatingFrameSunEarth(y_clean);
yrot_clean_halo = RotatingFrameSunEarth(y_clean_halo);

%%{
figure(3); %inertial frame
hold on
Moon = plot3(y_clean(:,1),y_clean(:,2),y_clean(:,3), 'k--' ,'DisplayName','Moon orbit');


[X,Y,Z] = sphere(20);
%L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
axis equal
%legend('show')
%legend(Moon)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
%}
%%

figure(4); %rotating frame
subplot(2, 3, 1);
[X,Y,Z] = sphere(20);
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'k+','DisplayName','Earth-Moon L2');
%Moon = plot3(rM,0,0,'go','DisplayName','Moon');
%plot3(y(:,1),y(:,2),y(:,3));  %plotting the Moon
Halo = plot3(yrot_clean_halo(:,7),yrot_clean_halo(:,8),yrot_clean_halo(:,9), 'k--' ,'DisplayName','Halo orbit around L1');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
Earth = plot3(0,0,0,'b.','MarkerSize',35,'DisplayName','Earth location')

%surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
Moon = plot3(rM,0,0,'g.','MarkerSize',25,'DisplayName','Moon location')
%legend([L1 Earth],{'Earth-Moon L1','Earth'})
%axis equal
%axis vis3d
%legend('show')
legend([Traj1 Traj2 Halo Insert_place L1 L2 Earth Moon ])
title('Zoomed in XY-view')

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

%%

figure(4); %rotating frame
subplot(2, 3, 2);
[X,Y,Z] = sphere(20);
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'k+','DisplayName','Earth-Moon L2');
%Moon = plot3(rM,0,0,'go','DisplayName','Moon');
%plot3(y(:,1),y(:,2),y(:,3));  %plotting the Moon
Halo = plot3(yrot_clean_halo(:,7),yrot_clean_halo(:,8),yrot_clean_halo(:,9), 'k--' ,'DisplayName','Halo orbit around L1');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
Earth = plot3(0,0,0,'b.','MarkerSize',35,'DisplayName','Earth location')

%surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
Moon = plot3(rM,0,0,'g.','MarkerSize',25,'DisplayName','Moon location')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Trajectories in XY-view')

%%

figure(4); %rotating frame
subplot(2, 3, 3);
[X,Y,Z] = sphere(20);
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'k+','DisplayName','Earth-Moon L2');
%Moon = plot3(rM,0,0,'go','DisplayName','Moon');
%plot3(y(:,1),y(:,2),y(:,3));  %plotting the Moon
Halo = plot3(yrot_clean_halo(:,7),yrot_clean_halo(:,8),yrot_clean_halo(:,9), 'k--' ,'DisplayName','Halo orbit around L1');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
Earth = plot3(0,0,0,'b.','MarkerSize',35,'DisplayName','Earth location')

%surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
Moon = plot3(rM,0,0,'g.','MarkerSize',25,'DisplayName','Moon location')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Trajectories in YZ-view')


%%
hold on
figure(4); %rotating frame
subplot(2, 3, 4);
[X,Y,Z] = sphere(20);
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'k+','DisplayName','Earth-Moon L2');
%Moon = plot3(rM,0,0,'go','DisplayName','Moon');
%plot3(y(:,1),y(:,2),y(:,3));  %plotting the Moon
Halo = plot3(yrot_clean_halo(:,7),yrot_clean_halo(:,8),yrot_clean_halo(:,9), 'k--' ,'DisplayName','Halo orbit around L1');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
Earth = plot3(0,0,0,'b.','MarkerSize',35,'DisplayName','Earth location')

%surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
Moon = plot3(rM,0,0,'g.','MarkerSize',25,'DisplayName','Moon location')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Zoomed in XZ-view')

%%
hold on
figure(4); %rotating frame
subplot(2, 3, 5);
[X,Y,Z] = sphere(20);
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'k+','DisplayName','Earth-Moon L2');
%Moon = plot3(rM,0,0,'go','DisplayName','Moon');
%plot3(y(:,1),y(:,2),y(:,3));  %plotting the Moon
Halo = plot3(yrot_clean_halo(:,7),yrot_clean_halo(:,8),yrot_clean_halo(:,9), 'k--' ,'DisplayName','Halo orbit around L1');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
Earth = plot3(0,0,0,'b.','MarkerSize',35,'DisplayName','Earth location')

%surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
Moon = plot3(rM,0,0,'g.','MarkerSize',25,'DisplayName','Moon location')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Trajectories in XZ-view')
 
%}
 
 %title(['Earth fixed rotating frame with manifold trajectories initiated by kick of ' num2str(dV_kick_p3) ' deltaV'])



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
y_axis_cross = 3.978726194287806e+08 - sqrt((y(7)^2)+(y(8)^2));
y_axis_cross = 4.0e+08 - sqrt((y(7)^2)+(y(8)^2));



value = [crashmoon; y_axis_cross];
isterminal = [1;1];
direction = [0;0];

end
