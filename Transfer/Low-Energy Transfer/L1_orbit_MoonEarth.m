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
r_halo_orbit = 38000000;
rot_speed_Moon = vM/rM;
v_L1 = rot_speed_Moon*(rM-rL1);

% Time
T_days = 12.1;  %T_days 12.1 for one orbit
T = 60*60*24*T_days; %s

% SUN
[X,Y,Z] = sphere(40);

%options of the integrator
options1 = odeset('RelTol', 1e-18);

[t1,y1] = ode113(@EarthMoonAcc,[0 T],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+896.5178218 0 ],options1); %initial orbit
yrot1 = RotatingFrameSunEarth(y1);

max_distance_X_direction = max(yrot1(:,7))-min(yrot1(:,7));
max_distance_Y_direction = max(yrot1(:,8))-min(yrot1(:,8));
max_distance_Z_direction = max(yrot1(:,9))-min(yrot1(:,9));


hold on
title(['L1 Earth-Moon halo orbit starting at z= ' num2str(r_halo_orbit) ' meter' ])

%plotting graph, bodies and Lagrane points
plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9),'DisplayName','halo orbit');
Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth,'DisplayName','Earth');
surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
L2 = plot3(rL2,0,0,'ko','DisplayName','Earth-Moon L2');

hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('show')
axis equal
axis vis3d