%Initial values
DOD_battery = 0.4;
P_req = 32; %Watt
N_cells = 30; %initial guess (sized in tool)
incl_constant = 50.2+6.68;

[T,T_eclipse] = Eclipse(1629e3,5);
T_sun = T-T_eclipse;
Temp_panel = ThermalData(incl_constant);
Irradiance = horzcat(ones(1,round(T_sun)),zeros(1,round(T_eclipse)));
inclination = cos(deg2rad(incl_constant))*ones(1,round(T));
E_batt = zeros(round(T)+1,1);
P_sol = zeros(round(T)+1,1);

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
        DT = Temp_panel(t)-temp_ref;
        P_sol(t) = (V_mp+DVDT*DT)*(I_mp+DIDT*DT)*N_cells*Irradiance(t)*inclination(t); %Power per cell*number of cells
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
        sizing = 0;
    elseif E_batt(end)-E_batt(1)<0
        N_cells = N_cells+1;
    elseif E_batt(end)-E_batt(1)>B_upperbound
        N_cells = N_cells-1;
    else
        error('Battery capacity not recognised')
    end
    
    
end

A_tot = N_cells*0.08*0.08;
P_inc = A_tot*1300;

BattCapacity = (max(E_batt)-min(E_batt))/(1-DOD_battery);
E_batt = E_batt + (BattCapacity-max(E_batt));
SOC_batt = E_batt./BattCapacity;
plot(0:T,E_batt)
% figure;
% plot(0:T,P_sol./P_inc)
notice_pwrsub = ['Solar Cells: ' num2str(N_cells) '. Solar panel area: '...
    num2str(A_tot) '. Battery capacity[Wh]: ' num2str(BattCapacity)];
