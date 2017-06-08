clear all;
close all;
clc;

Mu_Moon = 4.902801e12; %m^3 s^-2
R_Moon = 1737400; %m mean Moon radius

Matrix = xlsread('C:\Users\Bart\Desktop\DSE\DSE_Moon_TBD\Deployment\DeploymentSimulation.xlsx');

InitCon = [];
for i = linspace(1,size(Matrix,1),size(Matrix,1))
    InitCon = [InitCon Matrix(i,:)];
end

%options of the integrator
options = odeset('RelTol', 1e-12);
T = 100000; %s

[t,y] = ode113(@DeploymentSimulationAcc, [0 T], InitCon, options);

[X,Y,Z] = sphere(40);

figure
hold on
surf(X*R_Moon,Y*R_Moon,Z*R_Moon)
for i = linspace(1,size(y,2)-5,(size(y,2)/6))
    if i ==1
        plot3(y(:,1),y(:,2),y(:,3),'b')
    else
        plot3(y(:,i),y(:,(i+1)),y(:,(i+2)),'g')
    end        
end
title('Simulation of deployment')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
axis equal
axis vis3d



