%kepler2spherical converter function

%parameters
R_M = 1737.1e3; %[m]
M_M = 7.342e22; %[kg]
G = 6.67408e-11;
GMm = G*M_M;
%% inputs
e = 0.5; %[-]         Eccentricity
i = deg2rad(30); %[deg]       Inclination
RAAN = deg2rad(0); %[deg]    Right argument of the ascending node
omega = deg2rad(0); %[deg] argument of pericenter
phase = 0; % initial phase starting point



h = 400e3; %[m]
R_S = R_M+h; %[m] minimum altitude (periapsis)
a = R_S/(1-e);

%kep2chart_2002

nu=(linspace(phase,2*pi+phase,100))';
r = a*(1-e^2)./(1+e.*cos(nu));

X = r.*(cos(RAAN).*cos(omega+nu)-sin(RAAN).*sin(omega+nu)*cos(i));
Y = r.*(sin(RAAN).*cos(omega+nu)+cos(RAAN).*sin(omega+nu)*cos(i));
Z = r.*(sin(i).*sin(omega+nu));
[p,q,r] = sphere(40);
p = p.*R_M; q = q.*R_M; r = r.*R_M;

surf(p,q,r);hold on;
plot3(X,Y,Z,'LineWidth',1.5);hold on;
scatter3(X(1),Y(1),Z(1),16);hold on;
axis vis3d
axis equal

