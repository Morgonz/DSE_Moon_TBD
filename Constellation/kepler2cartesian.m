%kepler2spherical converter function

%% inputs
a = 0; %[m]        Semi-major axis
e = 0; %[-]         Eccentricity
i = 0; %[deg]       Inclination
RAAN = deg2rad(0); %[deg]    Right argument of the ascending node
%phasing is done automatically with uniform satellite spacing

%parameters
R_M = 1737.1e3; %[m]
H_S = 400e3; %[m]
R_S = R_M+H_S; %[m]

M_M = 7.342e22; %[kg]
G = 6.67408e-11;
mu_M = G*M_M;

%start position
x = R_S*cos(RAAN);
y = R_S*sin(RAAN);
z = 0;

d = [x;y;z];

T1 = [1 0 0;...
      0 cos(i) sin(i);...
      0 sin(i) cos(i)];
T2
phase = 0;
for a=(0+phase:pi/10:(2*pi+phase)
    z = [z 
    %calculate position
end

%% Outputs
R(t)
theta(t)
phi(t)
