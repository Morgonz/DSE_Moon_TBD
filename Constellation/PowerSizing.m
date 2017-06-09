%function input
h = 1629e3; % [m]
P_req = 32; % [W]
P_inc = 1300; % [W/m^2] incoming power
lifetime = 5; % Years
degradation = 0.03; % 0.0525 from SMAD, 0.03 from badescu2012moon
absorp = 0.9; % Solar absorptance

%Moon parameters
M_moon = 7.342e22; % [kg]
G = 6.67408e-11;
Rm = 1738.1e3; % [m]
% d_SE = 149597871e3; % [m]

GMm = G*M_moon;
Rs = Rm + h; % [m]

%eclipse time
T = 2*pi*sqrt(Rs^3/GMm);
eclipse_angle = 2*(asin(Rm/Rs));
T_eclipse = eclipse_angle*sqrt(Rs^3/GMm);
%T_eclipse = 0; %for halo sats

cycles = (lifetime*365*24*3600)/T;


%% Power calculations
%efficiency parameters
DOD_bat = 0.60; % Depth of discharge (Lithium Ion, 9e3 cycles)
%1000 to 25 degrees, 98% energy efficiency
e_sol = 0.26; %solar panel efficiency at BOL(triple junction GaAs/GaAs (improved MJ))
% Cell = 4x6 cm, -80 to 100 degrees, 0.21%/K (badescu2012moon)
L_D = (1-degradation)^lifetime; % degradation factor
L_P = cos(deg2rad(50.2+6.68)); % pointing factor: inclination + Moon tilt worst case


e_night = 0.65; %battery charging * discharging efficiency
e_day = 0.85; %battery charging

%Power calc (simulink integration?)
P_gen = ((P_req*(T-T_eclipse))/e_day+P_req*T_eclipse/e_day)/(T-T_eclipse); %SMAD
A_sol = P_gen/(P_inc*e_sol);% energy receival required
E_bat = P_req * T_eclipse / e_night / DOD_bat; % battery size requirement
C_bat = E_bat/3600; % [Wh]

% Area increase due to degradation and pointing offset
A_sol = A_sol / (L_D*L_P*absorp);



M_batt = C_bat/125; % from SMAD LI_ION battery
M_panel = A_sol*2.8;% from fake SMAD GaAs

M = M_batt + M_panel;

notice_pwr = ['A_panel[m^2]= ' num2str(A_sol) ' C_bat [Wh]= ' num2str(C_bat) ' Cycles: ' num2str(cycles) ' Mass [kg]: ' num2str(M)];
disp(notice_pwr)