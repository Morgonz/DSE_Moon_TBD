function [X,Y,Z] = keplerplot(in)
% IN = e,i,RAAN,omega,phase,h
e = in(1); i = in(2); RAAN = in(3); omega = in(4); phase = in(5); h = in(6);
%keplerplot: Translate kepler input into 3d plot
%   Detailed explanation goes here
%parameters
R_M = 1737.1e3; %[m]

R_S = R_M+h; %[m] minimum altitude (periapsis)
a = R_S/(1-e);

% e = 0.5; %[-]         Eccentricity
% i = deg2rad(30); %[deg]       Inclination
% RAAN = deg2rad(0); %[deg]    Right argument of the ascending node
% omega = deg2rad(0); %[deg] argument of pericenter
% phase = 0; % initial phase starting point
% h_per = 400e3; %[m]

%kep2chart_2002

theta=(linspace(phase,2*pi+phase,100))';
r = a*(1-e^2)./(1+e.*cos(theta));

X = r.*(cos(RAAN).*cos(omega+theta)-sin(RAAN).*sin(omega+theta)*cos(i));
Y = r.*(sin(RAAN).*cos(omega+theta)+cos(RAAN).*sin(omega+theta)*cos(i));
Z = r.*(sin(i).*sin(omega+theta));
% [p,q,r] = sphere(40);
% p = p.*R_M; q = q.*R_M; r = r.*R_M;
% 
% surf(p,q,r);hold on;
% plot3(X,Y,Z,'LineWidth',1.5);hold on;
% scatter3(X(1),Y(1),Z(1),16);hold on;
% axis vis3d
% axis equal

end

