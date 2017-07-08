clear all;
close all;
clc;

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Moon = 4.902801e12; %m^3 s^-2
R_Earth = 6371000; %m mean radius of the Earth
R_Moon = 1737400; %m mean Moon radius

%Circular Orbit initial conditions
h = 500000; %m
r1 = R_Earth + h; %m
V0 = sqrt(Mu_Earth/(r1)); %m/s

%Moon initial conditions
rM = 385000600; %m
vM = sqrt((Mu_Earth+Mu_Moon)/rM); %m/s speed of rotation of the Moon

%Earth-Moon L1, seen from the Moon
rHill = rM*(Mu_Moon/(3*Mu_Earth))^(1/3);

%Satellite intital conditions
rL1 = rM-rHill;
rL2 = rM+rHill;
r_halo_orbit = 44000000;
rot_speed_Moon = vM/rM;
v_L1 = rot_speed_Moon*(rL1);

% Time
T_days = 4*12.25-0.043362;  %T_days 12.25 for one orbit
T = 60*60*24*T_days; %s

% SUN
[X,Y,Z] = sphere(40);

%options of the integrator
options1 = odeset('RelTol', 2.22045e-14);

[t1,y1] = ode113(@EarthMoonAcc,[0 T],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1); %initial orbit
yrot1 = RotatingFrameSunEarth(y1);

%mirrored orbit
[t11,y11] = ode113(@EarthMoonAcc,[0 T],[rM 0 0 0 vM 0  rL1 0 -r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1); %initial orbit
yrot11 = RotatingFrameSunEarth(y11);


%[t1,y1] = ode113(@EarthMoonAcc,[0 T],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+234.91105658342964091 0 ],options1); %initial orbit
T_manouvre = 2*3600
diff = yrot1(end,:)-yrot1(1,:)
diff_pos = sqrt(diff(7)^2+diff(8)^2+diff(9)^2)
dV_needed = diff_pos/(T_manouvre*10)
dv_x = (y1(end,7)/diff_pos)*dV_needed;
dv_y = (y1(end,8)/diff_pos)*dV_needed;
dv_z = (y1(end,9)/diff_pos)*dV_needed;

[t2,y2] = ode113(@EarthMoonAcc,[0 0.32],[y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+dv_x y1(end,11)+dv_y y1(end,12)+dv_z],options1); %initial orbit
yrot2 = RotatingFrameSunEarth(y2);

dv_x_2 = yrot1(1,7)-yrot2(end,7)
dv_y_2 = yrot1(1,8)-yrot2(end,8)
dv_z_2 = yrot1(1,9)-yrot2(end,9)
V_diff_2 = sqrt(dv_x_2^2+dv_y_2^2+dv_z_2^2)

[t3,y3] = ode113(@EarthMoonAcc,[0 T/4],[y2(end,1) y2(end,2) y2(end,3) y2(end,4) y2(end,5) y2(end,6) y2(end,7) y2(end,8) y2(end,9) y2(end,10)+dv_x_2 y2(end,11)+dv_y_2 y2(end,12)+dv_z_2],options1); %initial orbit



%{
figure(1)
hold on
title('L1 Earth-Moon halo orbit in Earth fixed inertial frame')
plot3(y1(:,7),y1(:,8),y1(:,9),'r','DisplayName','halo orbit');
plot3(y1(:,1),y1(:,2),y1(:,3),'DisplayName','Moon orbit');
Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
Begin_Moon = plot3(y1(1,1),y1(1,2),y1(1,3),'ko','DisplayName','Start Moon');
Begin_Sat = plot3(y1(1,7),y1(1,8),y1(1,9),'ro','DisplayName','Start Satellite');
End_Moon = plot3(y1(end,1),y1(end,2),y1(end,3),'k*','DisplayName','Start Moon');
End_Sat = plot3(y1(end,7),y1(end,8),y1(end,9),'r*','DisplayName','Start Satellite');

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('show')
axis equal
axis vis3d

hold off
%}

%yrot1 = RotatingFrameSunEarth(y1);
%yrot2 = RotatingFrameSunEarth(y2);
yrot3 = RotatingFrameSunEarth(y3);

max_distance_X_direction = max(yrot1(:,7))-min(yrot1(:,7));
max_distance_Y_direction = max(yrot1(:,8))-min(yrot1(:,8));
max_distance_Z_direction = max(yrot1(:,9))-min(yrot1(:,9));


%figure(2)
hold on
figure(1)
subplot(2, 2, 1);       % add first plot in 2 x 2 grid
hold on
plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9),'DisplayName','halo orbit');
plot3(yrot11(:,7),yrot11(:,8),yrot11(:,9),'DisplayName','mirrored halo orbit');         % line plot
Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'ko','DisplayName','Earth-Moon L2');

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('show')
axis equal
axis vis3d

title('Halo orbits XY-view')


hold on
subplot(2,2,2)       % add second plot in 2 x 2 grid
hold on
plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9),'DisplayName','halo orbit');
plot3(yrot11(:,7),yrot11(:,8),yrot11(:,9),'DisplayName','mirrored halo orbit');        % scatter plot
Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'ko','DisplayName','Earth-Moon L2');

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
%legend('show')
axis equal
axis vis3d

title('Halo orbits YZ-view')


hold on
subplot(2,2,3)       % add third plot in 2 x 2 grid
hold on
plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9),'DisplayName','halo orbit');
plot3(yrot11(:,7),yrot11(:,8),yrot11(:,9),'DisplayName','mirrored halo orbit');          % stem plot
Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'ko','DisplayName','Earth-Moon L2');

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
%legend('show')
axis equal
axis vis3d


title('Halo orbits XZ-view')
%tightfig;
%}









%hold on
%title(['L1 Earth-Moon halo orbit starting at z= ' num2str(r_halo_orbit) ' meter' ])
%title('L1 Earth-Moon halo orbits')

%plotting graph, bodies and Lagrane points
%plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9),'DisplayName','halo orbit');
%plot3(yrot11(:,7),yrot11(:,8),yrot11(:,9),'DisplayName','mirrored halo orbit');
%plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9),'DisplayName','halo orbit maintenance');
%plot3(yrot3(:,7),yrot3(:,8),yrot3(:,9),'DisplayName','halo orbit after maintenance');
%End_halo = plot3(yrot1(end,7),yrot1(end,8),yrot1(end,9),'bo','DisplayName','end of halo orbit');
%Begin_halo = plot3(yrot1(1,7),yrot1(1,8),yrot1(1,9),'bs','DisplayName','begin of halo orbit');
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
%surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
%L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
%L2 = plot3(rL2,0,0,'ko','DisplayName','Earth-Moon L2');
%{
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('show')
axis equal
axis vis3d
%}
%hold off


%diff = yrot1(end,:)-yrot1(1,:)
%diff_pos = sqrt(diff(7)^2+diff(8)^2+diff(9)^2)
%diff_vel = sqrt(diff(10)^2+diff(11)^2+diff(12)^2)









