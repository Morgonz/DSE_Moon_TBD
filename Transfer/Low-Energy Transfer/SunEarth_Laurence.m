clear all;
close all;
clc;

%defining Mu from Nasa site (all 4 from: https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
Mu_Earth = 0.39860E+15
Mu_Sun = 132712E+15
R_Earth = 6371000; %m
R_Sun = 695700000; %m  

%Earth initial conditions
rE = 149.60E+9; %m ditance Earth to Sun from: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
vE = sqrt((Mu_Earth+Mu_Sun)/rE); %rotational speed of the Earth

%Satellite initial conditions aroun the Earth
h = 500000; %m
vS = sqrt(Mu_Earth/(R_Earth + h)); %m/s
rS = R_Earth + h; %m

%Sun-Earth L2 seen from the Earth
rL = rE*(Mu_Earth/(3*Mu_Sun))^(1/3);

%options of the integrator
options = odeset('RelTol', 1e-20);
T = 30000000; %s

[t,y] = ode113(@SunEarthAcc, [0 T], [rE 0 0 0 vE 0 rE+rL 0 0 0 (vE/rE)*(rE+rL) 0],options);

% SUN
[X,Y,Z] = sphere(40);

%plotting the figure
figure
hold on
axis equal

%Sun reference frame
plot3(y(:,1),y(:,2),y(:,3),'b','DisplayName','Earth trajectory')
plot3(y(:,7),y(:,8),y(:,9),'g','DisplayName','Satellite trajectory')
surf(X*R_Sun,Y*R_Sun, Z*R_Sun,'DisplayName','Sun position')

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