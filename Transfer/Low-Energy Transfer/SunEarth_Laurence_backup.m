clear all;
close all;
clc;

%defining Mu from Nasa site (all 4 from: https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
%Mu_Earth = 0.39860E+15;
%Mu_Sun = 132712E+15;
G =  6.67408E-11;
Mass_Earth = 5.9724E+24;
Mass_Sun = 1988500E+24;
Mu_Earth = Mass_Earth*G;
Mu_Sun = Mass_Sun*G;
%R_Earth = 6371000; %m
%R_Sun = 695700000; %m  

%Earth initial conditions
rE = 149.60E+9; %m ditance Earth to Sun from: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
vE = sqrt((Mu_Sun+Mu_Earth)/rE); %rotational speed of the Earth

r_baricentre = rE*(Mass_Earth/(Mass_Sun+Mass_Earth))

%Sun-Earth L2 seen from the Earth
%rL_Hill = rE*(Mu_Earth/(3*Mu_Sun))^(1/3);

% distance Sun-Earth L2 seen from the Sun
%dL2 = 1.511E+11;   %from https://en.wikipedia.org/wiki/Lagrangian_point
%vL2_rot = (vE/rE)*dL2;

r_L2_Earth = 1.50155E+9;  %calculated with wolframalpha
r_L2_Sun = rE+r_L2_Earth
v_rot_L2 = (vE/(rE-r_baricentre))*(r_L2_Sun-r_baricentre);


%options of the integrator
options = odeset('RelTol', 1e-15);
T = 4*25000000; %s

[t,y] = ode113(@SunEarthAcc, [0 T], [rE 0 0 0 vE 0 r_L2_Sun 0 0 0 v_rot_L2+5000 0],options);

%[t,y] = ode113(@SunEarthAcc, [0 T], [rE 0 0 0 vE 0 dL2 0 0 0 vL2_rot_L2 0],options);




yrot = RotatingFrameSunEarth(y);

% SUN
%[X,Y,Z] = sphere(40);

%plotting the figure
figure
hold on
axis equal

%Sun reference fixed frame
%plot3(y(:,1),y(:,2),y(:,3),'g','DisplayName','Earth trajectory')
%plot3(y(:,7),y(:,8),y(:,9),'b','DisplayName','Satellite trajectory')

%Sun reference rotating frame
plot3(yrot(:,1),yrot(:,2),yrot(:,3),'g*','DisplayName','Earth trajectory')
plot3(yrot(:,7),yrot(:,8),yrot(:,9),'b','DisplayName','Satellite trajectory')

%surf(X*R_Sun,Y*R_Sun, Z*R_Sun,'DisplayName','Sun position')

%Earth reference frame
%plot3(y(:,7)-y(:,1),y(:,8)-y(:,2),y(:,9)-y(:,3),'g')
%plot3(rL,0,0,'c*','DisplayName','La grange position')


%Earth-Sun L2 frame
%plot3(y(:,7)-(rE+rL),y(:,8)-(rE+rL),y(:,9)-(rE+rL),'g','DisplayName','Satellite trajectory')


 

hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('show')

%axis vis3d