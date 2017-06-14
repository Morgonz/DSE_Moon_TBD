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
vM = sqrt((Mu_Earth+Mu_Moon)/rM); %m/s speed of the Moon

%Earth-Moon L1, seen from the Moon
rHill = rM*(Mu_Moon/(3*Mu_Earth))^(1/3);

%Satellite intital conditions
%rL1 = rM-rHill;
rL2 = rM+rHill;
r_halo2_orbit = 42000000;
rot_speed_Moon = vM/rM; %rad/s
%v_L1 = rot_speed_Moon*rL1;
v_L2 = rot_speed_Moon*(rL2);
v_L2_orbit = rot_speed_Moon*(rL2+r_halo2_orbit);

%options of the integrator
options1 = odeset('RelTol', 1e-18);
options2 = odeset('RelTol', 1e-18, 'Events', @toocloseEarth);

%define delta v of kick adn total time
dV_kick = -100; %in delta V
T_days = 50;
T = 60*60*24*T_days; %s

%integrator

%
for t_manouvre = linspace(0.1,15.1,151)
    t = 60*60*24*t_manouvre; %s
    
    [t1,y1] = ode113(@EarthMoonAcc,[0 -t],[rM 0 0 0 vM 0  rL2+2000000 0 r_halo2_orbit 0 v_L2_orbit-240.65524714247 0],options1); %initial orbit
    yrot1 = RotatingFrameSunEarth(y1);
    angle = atan2(yrot1(end,8),(yrot1(end,7)-(rL2)));
    dv_x = dV_kick*sin(angle);
    dv_y = dV_kick*-cos(angle);
    
    [t2,y2] = ode113(@EarthMoonAcc, [-t -T], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7)+dV_kick*sin(angle) y1(end,8)-dV_kick*cos(angle) y1(end,9) y1(end,10) y1(end,11) y1(end,12)],options1);
    yrot2 = RotatingFrameSunEarth(y2);
    
    %plotting of the trajectories
    
    if yrot2(end,7)>rL2*1.2
        hold on
        %%% in Eearth fixed rotating frame
        %plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9));
        %plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9));

        %%% in Earth fixed inertial frame
        plot3(y1(:,1),y1(:,2),y1(:,3),'DisplayName','Moon trajectory'); %plot moon trajectory
        plot3(y1(:,7),y1(:,8),y1(:,9),'DisplayName','Sat trajectory'); %plot satellite trajectory
        plot3(y2(:,1),y2(:,2),y2(:,3),'DisplayName','Moon trajectory'); %plot moon trajectory
        plot3(y2(:,7),y2(:,8),y2(:,9),'DisplayName','Sat trajectory'); %plot satellite trajectory
    end
   
end
%}

theta = (0/180)*pi

hold on
% In the Lagrange point stable, not in orbit though
% Model the halo orbit clockwise direction
%[t1,y1] = ode113(@EarthMoonAcc,[0 T],[rM 0 0 0 vM 0  rL2+2000000 0 r_halo2_orbit 0 v_L2_orbit-240.65524714247 0],options1); %initial orbit

%34.153261616035

%yrot1 = RotatingFrameSunEarth(y1);
%plot3(yrot1(:,7),yrot1(:,8),yrot1(:,9));

%%%intertial plotting
%plot3(y1(:,1),y1(:,2),y1(:,3),'DisplayName','Moon trajectory'); %plot moon trajectory
%plot3(y1(:,7),y1(:,8),y1(:,9),'DisplayName','Sat trajectory'); %plot satellite trajectory

[X,Y,Z] = sphere(20);
L2 = plot3(rL2,0,0,'k*','DisplayName','Earth-Moon L2');
%plot3(y(:,1),y(:,2),y(:,3));  %plotting the Moon
%Earth = surf(X*R_Earth,Y*R_Earth,Z*R_Earth);
surf(X*R_Moon+rM,Y*R_Moon,Z*R_Moon,'DisplayName','Moon');


hold off
%legend([L2 Earth],{'Earth-Moon L2','Earth'})
%legend('show')
axis equal
axis vis3d
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')


%set(gca,'Ydir','reverse')
%set(gca,'Xdir','reverse')

title('Rotating frame with Earth-Moon L2 point')


