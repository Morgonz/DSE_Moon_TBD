%Eclipsetime
%finput:h, P_req

%% Eclipse calc
%Moon parameters
M_moon = 7.342e22; % [kg]
M_earth = 5.97219e24; % [kg]
G = 6.67408e-11;
R_M = 1737.1e3; % [m]
R_E = 6378.136e3; % [m]
d_EM = 384400e3; % [m]

GMe = G*M_earth;
GMm = G*M_moon;

h = 110e3; % [m]
R_S = R_M + h; % [m]

%Full orbit time
T = 2*pi*sqrt(R_S^3/GMm);

%simple case: Only Moon eclipse
e_angle = 2*(asin(R_M/R_S));
T_e = e_angle*sqrt(R_S^3/GMm);

SER_simple = T_e/T

%Worst case: Moon eclipse into Earth eclipse with one additional Moon eclipse
ea_angle = 2*atan(R_M/d_EM);

T_ea = ea_angle*sqrt(d_EM^3/GMe);

MAX_eclipse = 2*T_e + T_ea;
MAX_sun = T-T_e-T_ea;

SER_wc = MAX_eclipse/MAX_sun

%% Power calculations
%calculate power dimensions with SER usage
DOD_bat = 0.8; % Depth of discharge
E_sol = 0.28; %solar panel efficiency
E_bat = 0.7*0.6; %battery charging * discharging efficiency

