%function input
h = 1629e3; % [m]
P_req = 83.8*1.3; % [W]
P_inc = 1300; % [W/m^2] incoming power
lifetime = 5; % Years
degradation = 0.00; % 0.0525 from SMAD, 0.03 from badescu2012moon
absorp = 0.90; % Solar absorptance

[T,T_eclipse,cycles] = Eclipse(h,lifetime);
T = 29.5*24*3600;
T_eclipse = 0;%ceil(0.15*24*3600); % For NO ECLIPSE

%% Power calculations
%efficiency parameters
DOD_bat = 0.90; % Depth of discharge (Lithium Ion, 9e3 cycles)
%10 to 25 degrees, 98% energy efficiency
e_sol = 0.29; %solar panel efficiency at BOL(triple junction GaAs/GaAs (improved MJ))
% Cell = 4x6 cm, -80 to 100 degrees, 0.21%/K (badescu2012moon)
L_D = (1-degradation)^lifetime; % degradation factor
L_P = cos(deg2rad(52+1.5)); % pointing factor: inclination + Moon tilt worst case


e_night = 0.65; %battery charging * discharging efficiency
e_day = 0.85; %battery charging

%Power calc
P_gen = ((P_req*(T-T_eclipse))/e_day+P_req*T_eclipse/e_day)/(T-T_eclipse); %SMAD


%% Instantaneous power calculations
impulse_pwr = 126; %Watt
impulse_t = 500;
ACC_DOD = 1; %acceptable DOD
inst_Ebat = (impulse_pwr*impulse_t/ACC_DOD)/3600
charge_t = 10*24*3600;
impulse_P_gen = inst_Ebat/charge_t/e_night*3600

% P_gen = P_gen + impulse_P_gen;
%% sec
A_sol = P_gen/(P_inc*e_sol);% energy receival required
E_bat = P_req * T_eclipse / e_night / DOD_bat; % battery size requirement

C_bat = E_bat/3600; % [Wh]
% if C_bat<inst_Cbat
%     C_bat = inst_Cbat;
% else
%     %empty
% end

% Area increase due to degradation and pointing offset
A_sol = A_sol / (L_D*L_P*absorp)*(8.5^2/60.36);



M_batt = C_bat/166.6; % from SMAD LI_ION battery
M_panel = A_sol*3.26; % from fake SMAD GaAs

M = (M_batt + M_panel)*1.25;

notice_pwr = ['A_panel[m^2]= ' num2str(A_sol) ' C_bat [Wh]= ' num2str(C_bat) ' Mass [kg]: ' num2str(M)];
disp(notice_pwr)