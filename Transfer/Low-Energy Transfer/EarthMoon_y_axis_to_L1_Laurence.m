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
rL1 = rM-rHill
r_halo_orbit = 38000000
rot_speed_Moon = vM/rM
v_L1 = rot_speed_Moon*(rM-rL1)


%options of the integrator
options1 = odeset('RelTol', 1e-18);
options2 = odeset('RelTol', 1e-18, 'Events', @toocloseEarth);

%define delta v of kick adn total time
dV_kick = -100 %in delta V
T_days = 50;
T = 60*60*24*T_days; %s

%integrator

for t_manouvre = linspace(0.1,12.1,121)
    t = 60*60*24*t_manouvre; %s
    
    [t1,y1] = ode113(@EarthMoonAcc,[0 -t],[rM 0 0 0 vM 0  rL1 0 r_halo_orbit 0 v_L1+896.5178218 0 ],options1); %initial orbit
    yrot1 = RotatingFrameSunEarth(y1);
    angle = atan2(yrot1(end,8),(yrot1(end,7)-(rL1)));
    dv_x = dV_kick*sin(angle);
    dv_y = dV_kick*-cos(angle);
    
    [t2,y2] = ode113(@EarthMoonAcc, [-t -T], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7)+dV_kick*sin(angle) y1(end,8)-dV_kick*cos(angle) y1(end,9) y1(end,10) y1(end,11) y1(end,12)],options1);
    yrot2 = RotatingFrameSunEarth(y2);
    
    if yrot2(end,7)>rL1
        hold on
        plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9));
        plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9));
    end
    
    

end

[X,Y,Z] = sphere(20);
L1 = plot3(rL1,0,0,'k*','DisplayName','Earth-Moon L1');
%plot3(y(:,1),y(:,2),y(:,3));  %plotting the Moon
Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth);
surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');
hold off
legend([L1 Earth],{'Earth-Moon L1','Earth'})
axis equal
%axis vis3d
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Rotating frame with Earth-Moon L1 point')


