%Inputs
T_sun = T - T_eclipse;
Temp_panel = ThermalData;
Irradiance = horzcat(ones(1,round(T_sun)),zeros(1,round(T_eclipse)));
inclination = cos(deg2rad(56))*ones(1,round(T));
E_batt = zeros(round(T)+1,1);
SOC_batt = zeros(round(T)+1,1);
P_sol = zeros(round(T)+1,1);
%Initial values
SOC_battery = 0.625;
P_req = 32; %Watt
BattCapacity = 70; %Wh
N_cells = 95;

temp_ref = 25; %K nominal operation
E_batt(1) = BattCapacity*SOC_battery;
SOC_batt(1) = SOC_battery;
%From azur space 80*80
%BOL at 1376 W/m^2
DVDT = -6.7e-3; %V/K at max power
DIDT = 0.24e-3; %A/K at max power
V_mp = 2411e-3; %Nominal voltage at mpp
I_mp = 1007e-3; %Nominal current at mpp
A_tot = N_cells*0.08*0.08;
P_inc = A_tot*1300;

for t=1:T
    DT = Temp_panel(t)-temp_ref;
    P_sol(t) = (V_mp+DVDT*DT)*(I_mp+DIDT*DT)*N_cells*Irradiance(t)*inclination(t); %Power per cell*number of cells
    P_batt = P_sol(t)-P_req;
    if P_batt<0 % P_req>P_sol Solar panel does not generate enough power
        E_batt(t+1) = E_batt(t)+P_batt/3600; %new capacity
    elseif P_batt>0 %Solar panel produces too much power
        E_batt(t+1) = E_batt(t)+P_batt/3600; %new capacity
    else %P_batt = 0, no active discharge
        E_batt(t+1) = E_batt(t)+P_batt/3600;
    end
    SOC_batt(t+1) = E_batt(t+1)/BattCapacity;
end
plot(0:T,SOC_batt)
figure;
plot(0:T,P_sol./P_inc)