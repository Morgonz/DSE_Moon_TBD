%Initial values
DOD_battery = 0.5; %capacity reduction to 70% at EOL
P_req = 45.40*1.3; %Watt + 0.3 MARGIN (smad)
N_cells = 72; %initial guess (sized in tool)
autosize = 1;
load solar_incidence
sat = 2;

%Initialise: network
incl_constant = 50.2+1.54;%53;
[T,T_eclipse] = Eclipse(1629e3,5);
T = round(T);%12.25*3600*24;
T_eclipse = round(T_eclipse);
T_sun = T-T_eclipse;
temp_ref = 28+273.15; %K nominal operation
Temp_panel = ThermalData(incl_constant,1629e3);
% Temp_panel = ones(size(Temp_panel))*(335.7);
Irradiance = horzcat(ones(1,round(T_sun)),zeros(1,round(T_eclipse)));
inclination = cos(deg2rad(incl_constant))*ones(1,round(T)); % cos(deg2rad(y));
mfac = 1;

% %Initialise: Relay
% T = 868;%*(12.25*3600*24/360);%(1:868).*(12.25*3600*24/360);
% temp_ref = 28+273.15; %K nominal operation
% Temp_panel = ones(1,868)*(335.7);
% Irradiance = ones(1,868);
% inclination = cos(deg2rad(y(sat,:)));
% mfac = 12.25*3600*24/360;


%% rest

E_batt = zeros(round(T)+1,1);
P_sol = zeros(round(T)+1,1);
DT = Temp_panel-temp_ref;
L_D = (1-0.03)^5;
absorbtance = 0.91;

B_upperbound = 2000;

E_batt(1) = 0;


%BOL at 1376 W/m^2
A_cell = 60.36;%cm^2
DVDT = -6.7e-3; %V/K at max power
DIDT = 0.24e-3; %A/K at max power
V_mp = 2.411; %Nominal voltage at mpp
I_mp = 1.007; %Nominal current at mpp
sizing = 1;


while sizing

    for t=1:T
        
        P_sol(t) = (V_mp+DVDT*DT(t))*(I_mp+DIDT*DT(t))*N_cells*Irradiance(t)*L_D*absorbtance*inclination(t); %Power per cell*number of cells
        P_batt = P_sol(t)-P_req;
        if P_batt<0 % P_req>P_sol Solar panel does not generate enough power
            E_batt(t+1) = E_batt(t)+P_batt/3600/0.65*mfac; %new capacity
        elseif P_batt>0 %Solar panel produces too much power
            E_batt(t+1) = E_batt(t)+P_batt/3600*0.85*mfac; %new capacity
        else %P_batt = 0, no active discharge
            E_batt(t+1) = E_batt(t)+P_batt/3600*mfac;
        end
    end
    
    if abs(E_batt(end)-E_batt(1)) <= B_upperbound && E_batt(end)-E_batt(1) >=0
        sizing = 0; %Solar panel sizing is within set bounds
    elseif E_batt(end)-E_batt(1)<0
        N_cells = N_cells+1; %Enlarge solar array
    elseif E_batt(end)-E_batt(1)>B_upperbound
        N_cells = N_cells-1; %Shrink solar array
    else
        error('Battery capacity not recognised') %?
    end
    if autosize == 0
        sizing = 0;
    end
    disp(N_cells)
end

A_tot = N_cells*A_cell/10000;
P_inc = A_tot*1300;

BattCapacity = (max(E_batt)-min(E_batt))/(DOD_battery);
E_batt = E_batt + (BattCapacity-max(E_batt));
SOC_batt = E_batt./BattCapacity;

%% Plot

% VerificationData = [0.1517502602 0.1357155525 0.1282785684 0.127196632 0.1271362953];
% VerificationPoint= [1 2000 5000 9000 14509];

figure;
plot(0:T,P_sol./P_inc)
figure;
yyaxis left
plot(0:T,P_sol,'LineWidth',1.5); hold on;
plot([0 T],[P_req P_req],'LineWidth',1.5);
ylabel('Power [W]')

yyaxis right
plot(0:T,SOC_batt,'LineWidth',1.5)
ylabel('Battery SOC [-]')
xlabel('Time [hr]')
xlabelloc = 0:3600:T;
xlim([0 T(end)])
xticks(xlabelloc)
xticklabels(0:length(xlabelloc))
legend('Solar panel power output','Power Requirement','Battery SOC','Location','South')



%% Sizing
%Battery
%Selected: Saft VL51ES
Batt_sp_energy = 170;
Batt_Volume = 2*pi*0.054*0.222;
V = 3.6;
Capacity = 51; % Ah
Mass_batt = 1.08; % kg5
Batt_energy = 180; % Wh
Batt_margin = Batt_energy/(BattCapacity/0.7);

%Solar cells
%Selected: Azur space 3g30c
mass_per_area = 1.3+0.59+1.37;
Mass_panel = mass_per_area * A_tot;
Mass_total = Mass_panel + Mass_batt;

notice_pwrsub = ['Solar Cells: ' num2str(N_cells) '. Solar panel area[m^2]: '...
    num2str(A_tot) '. Battery charge[Wh]: ' num2str(BattCapacity) '. System mass[kg]: ' num2str(Mass_total)];
disp(notice_pwrsub)
