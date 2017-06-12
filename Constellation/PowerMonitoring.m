%Initial values
DOD_battery = 0.4;
P_req = 32; %Watt
N_cells = 30; %initial guess (sized in tool)
incl_constant = 50.2+6.68;

[T,T_eclipse] = Eclipse(1629e3,5);
T = round(T);
T_eclipse = round(T_eclipse);
T_sun = T-T_eclipse;
Temp_panel = ThermalData(incl_constant,1629e3);
% Temp_panel = ones(size(ThermalData(incl_constant,1629e3)))*120;
Irradiance = horzcat(ones(1,round(T_sun)),zeros(1,round(T_eclipse)));
inclination = cos(deg2rad(incl_constant))*ones(1,round(T));
E_batt = zeros(round(T)+1,1);
P_sol = zeros(round(T)+1,1);
DT = Temp_panel-temp_ref;
L_D = (1-0.03)^5;
absorbtance = 0.9;

B_upperbound = 6;
temp_ref = 28+273.15; %K nominal operation
E_batt(1) = 0;

%From azur space 80*80
%BOL at 1376 W/m^2
DVDT = -6.7e-3; %V/K at max power
DIDT = 0.24e-3; %A/K at max power
V_mp = 2411e-3; %Nominal voltage at mpp
I_mp = 1007e-3; %Nominal current at mpp
sizing = 1;


while sizing
    disp(N_cells)
    for t=1:T
        P_sol(t) = (V_mp+DVDT*DT(t))*(I_mp+DIDT*DT(t))*N_cells*Irradiance(t)*L_D*absorbtance*inclination(t); %Power per cell*number of cells
        P_batt = P_sol(t)-P_req;
        if P_batt<0 % P_req>P_sol Solar panel does not generate enough power
            E_batt(t+1) = E_batt(t)+P_batt/3600/0.65; %new capacity
        elseif P_batt>0 %Solar panel produces too much power
            E_batt(t+1) = E_batt(t)+P_batt/3600*0.85; %new capacity
        else %P_batt = 0, no active discharge
            E_batt(t+1) = E_batt(t)+P_batt/3600;
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
end

A_tot = N_cells*0.08*0.08;
P_inc = A_tot*1300;

BattCapacity = (max(E_batt)-min(E_batt))/(DOD_battery);
E_batt = E_batt + (BattCapacity-max(E_batt));
SOC_batt = E_batt./BattCapacity;
%% Plot
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
legend('Solar panel','Requirement','Battery SOC','Location','NorthWest')

notice_pwrsub = ['Solar Cells: ' num2str(N_cells) '. Solar panel area: '...
    num2str(A_tot) '. Battery charge[Wh]: ' num2str(BattCapacity)];
disp(notice_pwrsub)