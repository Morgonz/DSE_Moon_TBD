%function input
h = 1629e3; % [m]
P_req = 32; % [W]
P_inc = 1300; % [W/m^2] incoming power
lifetime = 5;
degradation = 0.0525; % SMAD

%Moon parameters
M_moon = 7.342e22; % [kg]
G = 6.67408e-11;
Rm = 1738.1e3; % [m]
d_SE = 149597871e3; % [m]

GMm = G*M_moon;
Rs = Rm + h; % [m]

%eclipse time
T = 2*pi*sqrt(Rs^3/GMm);
cycles = (lifetime*365*24*3600)/T;
e_angle = 2*(asin(Rm/Rs));
T_eclipse = e_angle*sqrt(Rs^3/GMm);
%T_eclipse = 0; %for halo sats

%% Power calculations
%efficiency parameters
DOD_bat = 1-0.40; % Depth of discharge
e_sol = 0.30; %solar panel efficiency
e_night = 0.65; %battery charging * discharging efficiency
e_day = 0.85; %battery charging

%Power calc (simulink integration?)
P_gen = ((P_req*(T-T_eclipse))/e_day+P_req*T_eclipse/e_day)/(T-T_eclipse); %SMAD
A_sol = P_gen/(P_inc*e_sol);% energy receival required
E_bat = P_req * T_eclipse / e_night / DOD_bat; % battery size requirement
C_bat = E_bat/3600; % [Wh]

%Efficiency
L_D = (1-degradation)^lifetime; % degradation factor
L_A = cos(deg2rad(50.2+6.68)); % pointing factor: inclination + Moon tilt worst case
A_sol = A_sol / (L_D*L_A);

M_batt = C_bat/125; % from SMAD LI_ION battery
M_panel = A_sol*2.8;% from fake SMAD GaAs
M = M_batt + M_panel;

notice_pwr = ['Scenario = ' scen ' A_panel= ' num2str(A_sol) ' C_bat [Wh]= ' num2str(C_bat) ' Cycles: ' num2str(cycles) ' Mass [kg]: ' num2str(M)];
disp(notice_pwr)
