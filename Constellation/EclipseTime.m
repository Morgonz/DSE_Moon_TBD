%Eclipsetime 2030-2045:
%Penumbral:13 times
%Partial eclipse: 8 times
%full eclipse:12 times
%total full eclipse time: 890 min
%total partial eclipse time: 3732 min
%https://web.archive.org/web/20070305183925/http://sunearth.gsfc.nasa.gov/eclipse/LEcat/LE2001-2100.html

%function input
h = 1629e3; % [m]
P_req = 32; % [W]
P_inc = 1300; % [W/m^2] incoming power
lifetime = 5*365*24*3600;

%% Eclipse calc
%Moon parameters
M_moon = 7.342e22; % [kg]
M_earth = 5.97219e24; % [kg]
G = 6.67408e-11;
Rm = 1738.1e3; % [m]
Re = 6378.136e3; % [m]
Rs = 695700e3; % [m]
d_EM = 384400e3; % [m]
d_SE = 149597871e3; % [m]

GMe = G*M_earth;
GMm = G*M_moon;


Rs = Rm + h; % [m]


%Full orbit time
T = 2*pi*sqrt(Rs^3/GMm);
cycles = (5*365*24*3600)/T;

%% simple case: Only Moon eclipse
e_angle = 2*(asin(Rm/Rs));
T_eclipse = e_angle*sqrt(Rs^3/GMm);
%T_eclipse = 0; %for halo sats
SER_simple = T_eclipse/T;

%% Penumbral case: Full penumbra pass (also viable for partial eclipses)
a_penumbral = atan((Rs+Re)/d_SE);
R_penumbral = (d_EM+d_SE)*tan(a_penumbral) - Rs;
p_angle = 2*atan(R_penumbral/d_EM);
T_penumbral = p_angle*sqrt(d_EM^3/GMe);

%eclipse -> [S] penumbral -> eclipse -> penumbral -> ecli[E]pse

PEN_eclipse = 3*T_eclipse;
PEN_sun = 0.5*(2*(T-T_eclipse));
PEN_time = PEN_eclipse + PEN_sun;

SER_PEN = PEN_eclipse/PEN_time;

freq_PEN = (8*T_penumbral/2 + 13*T_penumbral/2)/lifetime;

%% Worst case: Moon eclipse into Earth eclipse with one additional Moon eclipse
ea_angle = 2*atan(Rm/d_EM);

T_ea = ea_angle*sqrt(d_EM^3/GMe);

%eclipse -> [B] solar ecl[E]ipse -> eclipse
MAX_eclipse = 2*T_eclipse + T_ea;
MAX_sun = T-T_eclipse-T_ea;
MAX_time = MAX_eclipse + MAX_sun;

SER_MAX = MAX_eclipse/MAX_time;

freq_MAX = (890*60)/lifetime;
%% Power calculations
%calculate power dimensions with SER usage
%Solar panel from https://www.isispace.nl/product/isis-cubesat-solar-panels/

%efficiency parameters
DOD_bat = 1-0.40; % Depth of discharge
e_sol = 0.30; %solar panel efficiency
e_night = 0.65; %battery charging * discharging efficiency
e_day = 0.85;

SER_USED = SER_simple;

%Power calc (simulink integration?)
E_req = P_req * T; %[J] energy per orbit
P_gen = ((P_req*(T-T_eclipse))/e_day+P_req*T_eclipse/e_day)/(T-T_eclipse);
%P_gen = ((P_req*(T-T_eclipse)+P_req*(T_eclipse)/e_night)/e_sol)/T;
A_sol = P_gen/(P_inc*e_sol);%((P_req/(1-SER_USED))/e_sol)/P_inc; % [W] energy receival required
E_bat = P_req * T_eclipse / e_night / DOD_bat; % battery size requirement
C_bat = E_bat/3600; % [Wh]

A_sol_wc = ((P_req/(1-SER_MAX))/e_sol)/P_inc;
E_bat_wc = P_req * SER_MAX*T / e_night / DOD_bat;

if SER_USED == SER_simple
    scen = 'simple.';
elseif SER_USED == SER_MAX
    scen = 'full eclipse.';
elseif SER_USED == SER_PEN
    scen = 'penumbral.';
else
    scen = 'unknown.';
end

degradation = 0.0525;
L_D = (1-degradation)^lifetime;
A_sol = A_sol / L_D/cos(deg2rad(50.2+6.68));

M_batt = C_bat/125; % from SMAD LI_ION battery
M_panel = A_sol*2.8;% from fake SMAD GaAs
M = M_batt + M_panel;

notice_pwr = ['Scenario = ' scen ' A_panel= ' num2str(A_sol) ' C_bat [Wh]= ' num2str(C_bat) ' DownFrac = ' num2str(freq_MAX+freq_PEN)];
disp(notice_pwr)