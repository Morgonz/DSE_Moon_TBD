%function input
h = 1629e3; % [m]
P_req = 104*1.3; % [W]
P_inc = 1300; % [W/m^2] incoming power
lifetime = 5; % Years
degradation = 0.03; % 0.0525 from SMAD, 0.03 from badescu2012moon
absorp = 0.91; % Solar absorptance

[T,T_eclipse,cycles] = Eclipse(h,lifetime);
%T_eclipse = 0; % For NO ECLIPSE

%% Power calculations
%efficiency parameters
DOD_bat = 0.50; % Depth of discharge (Lithium Ion, 9e3 cycles)
%1000 to 25 degrees, 98% energy efficiency
e_sol = 0.29; %solar panel efficiency at BOL(triple junction GaAs/GaAs (improved MJ))
% Cell = 4x6 cm, -80 to 100 degrees, 0.21%/K (badescu2012moon)
L_D = (1-degradation)^lifetime; % degradation factor
L_P = cos(deg2rad(50.2+1.5)); % pointing factor: inclination + Moon tilt worst case


e_night = 0.65; %battery charging * discharging efficiency
e_day = 0.85; %battery charging

%Power calc
P_gen = ((P_req*(T-T_eclipse))/e_day+P_req*T_eclipse/e_day)/(T-T_eclipse); %SMAD


%Instantaneous power calculations
impulse_pwr = 0; %Watt
impulse_t = 180;
ACC_DOD = 0.8; %acceptable DOD
inst_Cbat = impulse_pwr*impulse_t/ACC_DOD;
charge_t = 20*60;
impulse_P_gen = inst_Cbat/charge_t/e_night

% P_gen = P_gen + impulse_P_gen;

A_sol = P_gen/(P_inc*e_sol);% energy receival required
E_bat = P_req * T_eclipse / e_night / DOD_bat; % battery size requirement

C_bat = E_bat/3600; % [Wh]
% if C_bat<inst_Cbat
%     C_bat = inst_Cbat;
% else
%     %empty
% end

% Area increase due to degradation and pointing offset
A_sol = A_sol / (L_D*L_P*absorp);



M_batt = C_bat/170; % from SMAD LI_ION battery
M_panel = A_sol*3.26;% from fake SMAD GaAs

M = M_batt + M_panel;

notice_pwr = ['A_panel[m^2]= ' num2str(A_sol) ' C_bat [Wh]= ' num2str(C_bat) ' Cycles: ' num2str(cycles) ' Mass [kg]: ' num2str(M)];
disp(notice_pwr)