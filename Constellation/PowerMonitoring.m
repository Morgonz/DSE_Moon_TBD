%Initial values
DOD_battery = 0.4; %capacity reduction to 70% at EOL
P_req = (28.5); %Watt + 0.3 MARGIN (smad)
N_cells = 48; %initial guess (sized in tool)
autosize = 0;   
sat = 6;

%Initialise: network
incl_constant = 50.2+1.54;%53;
[T,T_eclipse] = Eclipse(1629e3,5);
T = round(T);%12.25*3600*24;
T_eclipse = round(T_eclipse);
T_sun = T-T_eclipse;
temp_ref = 28+273.15; %K nominal operation
Temp_panel = ThermalData(incl_constant,1629e3);
% Temp_panel = ones(size(Temp_panel)).*320;
Irradiance = horzcat(ones(1,round(T_sun)),zeros(1,round(T_eclipse)));
% Irradiance = ones(1,round(T));
incidence = cos(deg2rad(incl_constant))*ones(1,round(T));
mfac = 1;
bs = 0;

% %Initialise: Relay
% load solar_incidence
% 
% T = 868;%*(12.25*3600*24/360);%(1:868).*(12.25*3600*24/360);
% temp_ref = 28+273.15; %K nominal operation
% Temp_panel = ones(1,868)*(320);
% Irradiance = ones(1,868);
% incidence = cos(deg2rad(y(sat,:)));
% mfac = 12.25*3600*24/360;
% 
% bs = 1;

%% rest

E_batt = zeros(round(T)+1,1);
SOC_batt = zeros(round(T)+1,1);
diff = zeros(round(T),1);
P_sol = zeros(round(T)+1,1);
DT = Temp_panel-temp_ref;
L_D = (1-0.03)^5;
absorbtance = 0.91;

B_upperbound = 500000*mfac;

E_batt(1) = 0;
SOC_batt(1) = 1;

%BOL at 1376 W/m^2
A_cell = 60.36;%cm^2
DVDT = -6.7e-3; %V/K at max power
DIDT = 0.24e-3; %A/K at max power
V_mp = 2.411; %Nominal voltage at mpp
I_mp = 1.007; %Nominal current at mpp
sizing = 1;


while sizing

    for t=1:T
        
        P_sol(t) = (V_mp+DVDT*DT(t))*(I_mp+DIDT*DT(t))*N_cells*Irradiance(t)*L_D*absorbtance*incidence(t); %Power per cell*number of cells
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

if bs == 0
    BattCapacity = (max(E_batt)-min(E_batt))/(DOD_battery);
    E_batt = E_batt + (BattCapacity-max(E_batt));
    SOC_batt = E_batt./BattCapacity;
elseif bs == 1
    BattCapacity = 540;
    for bc=1:length(E_batt)-1
        diff(bc) = (E_batt(bc+1)-E_batt(bc))/BattCapacity;
        if diff(bc)>0 && SOC_batt(bc)>=1
            SOC_batt(bc+1) = 1;
        elseif diff(bc)>0 && (SOC_batt(bc)+diff(bc))<=1
            SOC_batt(bc+1) = SOC_batt(bc)+diff(bc);
        elseif diff(bc)>0 && (SOC_batt(bc)+diff(bc))>1
            SOC_batt(bc+1) = 1;
        elseif diff(bc)<=0
            SOC_batt(bc+1) = SOC_batt(bc)+diff(bc);
        else
            error('erreur')
        end
    end
else
    error('Battery sizing mode not recognised')
end
figure;
plot(SOC_batt); hold on;
plot(diff)


%% Plot

% VerificationData = [0.1517502602 0.1357155525 0.1282785684 0.127196632 0.1271362953];
% VerificationPoint= [1 2000 5000 9000 14509];
% figure;
% plot(1:T,y(1,:),1:T,y(3,:),1:T,y(4,:)',1:T,y(5,:)',1:T,y(6,:)');hold on;
% scatter(1:T,y(2,:),25,'filled')
% xlim([0 T])
% xlabel('Time [days]')
% ylabel('Solar angle of incidence')
% xticks(0:70:T)%0:3600:T
% xticklabels(0:floor(T/70))%T/3600
% set(gca,'FontSize',12)
% figure;
% yyaxis left
% plot(0:T,P_sol./P_inc,'LineWidth',1.5);
% ylabel('Solar panel efficiency [-]','FontSize',16)
% yyaxis right
% plot(1:T,incidence,'LineWidth',1.5);
% ylabel('Solar incidence factor [-]','FontSize',16)
% xlabel('Time [days]','FontSize',16)
% xlim([0 T])
% xticks(0:70:T)%0:3600:T
% xticklabels(0:floor(T/70))%T/3600
% set(gca,'FontSize',12)

figure;
yyaxis left
plot(0:T,P_sol./P_inc,'LineWidth',1.5);
ylabel('Solar panel efficiency [-]','FontSize',16)
yyaxis right
plot(1:T,Temp_panel-273.15,'LineWidth',1.5);
ylabel('Solar panel temperature [^\circC]','FontSize',16)
xlabel('Time [hr]','FontSize',16)
xlim([0 T])
xticks(0:3600:T)%0:3600:T
xticklabels(0:floor(T/3600))%T/3600
set(gca,'FontSize',12)



figure;
yyaxis left
plot(0:T,P_sol,'LineWidth',1.5); hold on;
plot([0 T],[P_req P_req],'LineWidth',1.5);
ylabel('Power [W]')

yyaxis right
plot(0:T,SOC_batt,'LineWidth',1.5)
ylabel('Battery SOC [-]')
xlabel('Time [h]')

xlim([0 T(end)])
xticks(0:3600:T)
xticklabels(0:floor(T/3600))
legend('Solar panel power output','Power Requirement','Battery state of charge','Location','South')
set(gca,'FontSize',12)


%% Sizing
% Battery
% Selected: Saft VES16

% amount_batt = 12;
% Batt_sp_energy = 155;
% Batt_Volume = pi*0.033^2*0.060;
% V = 3.6;
% Capacity = 4.5; % Ah
% Mass_batt = 0.155*amount_batt; % kg5
% Batt_energy = 16*amount_batt; % Wh
% Batt_margin = Batt_energy/(BattCapacity/0.7);
% disp(Batt_margin)
Mass_batt = BattCapacity/155;

%Solar cells
%Selected: Azur space 3g30c
mass_per_area = 1.3+0.59+1.37;
Mass_panel = mass_per_area * N_cells*0.085^2;
Mass_total = (Mass_panel + Mass_batt)*1.15;

mass_hinge = 0.01;

notice_pwrsub = ['Solar Cells: ' num2str(N_cells) '. A_panel[m^2]: '...
    num2str(A_tot) '. E_batt[Wh]: ' num2str(BattCapacity) ...
    '. System mass[kg]: ' num2str(Mass_total)];
disp(notice_pwrsub)