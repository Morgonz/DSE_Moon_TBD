function [ Temps_out ] = ThermalData(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Thermal controll program
 
 %Input parameters
 panel_incl = 0; % inclination (0 -> 50.2+6.68)
 T_i=288; %K initial temperature
 h=1629; %Orbit height moon (km)
 Period= round(1.753722994309795e+04); %Orbital period (s)
 Eclipse_T= round(3.027876566028392e+03); %Eclipse time (s)
 
 %% Selection of over how many orbits
 n=3; %Ammount of orbits
 
 
 %% L1-L2 lagrange points
 
 L1=148100000; %km
 L2=151100000; %km
 Sun_Earth=149600000; %km
 DL1_E=abs(Sun_Earth-L1); %km distance L1 to earth 
 DL2_E=abs(Sun_Earth-L2); %km Distance L2 to earth
 
 
 %% fixed parameters
 Roe=384400; %Radius orbit earth (km)
 Rm=1737;%radius moon (km)
 Re=6371;%radius earth (km)
 Rom=Rm+h; %Radius orbit moon
 Ae=0.306;%albedo earth
 Am=0.11;%albedo moon
 SF=1361;%W/m^2 solar constant
 L1SF=1361*(Sun_Earth^2/L1^2);%W/m^2 solar constant
 L2SF=1361*(Sun_Earth^2/L2^2);%W/m^2 solar constant
 EF=236*((Re/Roe))^2; %W/m^2 earth radiation
 MF=303*((Rm/Rom))^2; %W/m^2 moon radiation
 sigma=5.67036713*10^-8; %W/(m^2*K^4)
 R_E=(SF*(Ae)*(Re*1000)^2)/((Roe*1000)^2); %W/m^2 reflected earth radiation 
 R_M=(SF*(Am)*(Rm*1000)^2)/((Rom*1000)^2); %W/m^2 reflected moon radiation 

 %% Absorptance and Emmitance table 
 AE=[0.21 0.04     % 1 Aluminum tape
     0.24 0.08     % 2 Polished aluminium
     0.14 0.05     % 3 Aluminized kapton (alu outside)
     0.95 0.90     % 4 Black paint(polyurethane)
     0.29 0.83     % 5 White paint(silicone)
     0.88 0.80     % 6 solar panel(GaAs)
     0.07 0.74     % 7 Optical solar reflector)
     0.40 0.63     % 8 Aluminized kapton(alu inside)
     0.14 0.90     % 9 White paint silicate 
     0.25 0.02     % 10 Goldized kapton (gold outside)  
     0.95 0.80     % 11 Black paint electrical conducting
     0.37 0.44     % 12 Silver paint, electrical conducting
     0.77 0.84];    % 13 Anodized aluminum 
 
 
 
 %% Node classification
 
 % outer skins. [1 A(m^2) |2 thickness(m) |3 width(m) |4 height(m) |5 Cp(J/(K*kg)) |6 absorbtance |7 emmisivity
 %|8 conductivity(W/(m*K)) |9 density(kg/m^3)]
 Xp=[0.06 0.002 0.3 0.2 910 AE(13,1) AE(13,2) 204 2830];
 Xm=[0.06 0.002 0.3 0.2 910 AE(13,1) AE(13,2) 204 2830];
 Yp=[0.06 0.002 0.3 0.2 910 AE(13,1) AE(13,2) 204 2830];
 Ym=[0.06 0.002 0.3 0.2 910 AE(13,1) AE(13,2) 204 2830];
 Zp=[0.04 0.002 0.2 0.2 910 AE(13,1) AE(13,2) 204 2830];
 Zm=[0.04 0.002 0.2 0.2 910 AE(13,1) AE(13,2) 204 2830];
 % Tanks
 
 %Batteries
 
 %Solar panels
 Sp=[0.16 0.004 0.8 0.2 910 AE(6,1) AE(6,2) 204 2830];
 Spb=AE(9,2);

 %Antenna
 
 Nodes=[Xp
        Xm
        Yp
        Ym
        Zp
        Zm];
    
 %Conductive paths (conductivety*area/(length*spesific heat*mass)
 Xp_Yp=1/((1/((Xp(2)*Xp(3)*Xp(8)/(Xp(4)/2))/(Xp(1)*Xp(2)*Xp(9)*Xp(5))))+(1/((Yp(2)*Yp(3)*Yp(8)/(Yp(4)/2))/(Yp(1)*Yp(2)*Yp(9)*Yp(5)))));   %Xp to Yp
 Xp_Ym=1/((1/((Xp(2)*Xp(3)*Xp(8)/(Xp(4)/2))/(Xp(1)*Xp(2)*Xp(9)*Xp(5))))+(1/((Ym(2)*Ym(3)*Ym(8)/(Ym(4)/2))/(Ym(1)*Ym(2)*Ym(9)*Ym(5)))));   %Xp to Ym
 Xp_Zp=1/((1/((Xp(2)*Xp(4)*Xp(8)/(Xp(3)/2))/(Xp(1)*Xp(2)*Xp(9)*Xp(5))))+(1/((Zp(2)*Zp(4)*Zp(8)/(Zp(3)/2))/(Zp(1)*Zp(2)*Zp(9)*Zp(5)))));   %Xp to Zp
 Xp_Zm=1/((1/((Xp(2)*Xp(4)*Xp(8)/(Xp(3)/2))/(Xp(1)*Xp(2)*Xp(9)*Xp(5))))+(1/((Zm(2)*Zm(4)*Zm(8)/(Zm(3)/2))/(Zm(1)*Zm(2)*Zm(9)*Zm(5)))));   %Xp to Zm
 
 Xm_Yp=1/((1/((Xm(2)*Xm(3)*Xm(8)/(Xm(4)/2))/(Xm(1)*Xm(2)*Xm(9)*Xm(5))))+(1/((Yp(2)*Yp(3)*Yp(8)/(Yp(4)/2))/(Yp(1)*Yp(2)*Yp(9)*Yp(5)))));   %Xm to Yp
 Xm_Ym=1/((1/((Xm(2)*Xm(3)*Xm(8)/(Xm(4)/2))/(Xm(1)*Xm(2)*Xm(9)*Xm(5))))+(1/((Ym(2)*Ym(3)*Ym(8)/(Ym(4)/2))/(Ym(1)*Ym(2)*Ym(9)*Ym(5)))));   %Xm to Ym
 Xm_Zp=1/((1/((Xm(2)*Xm(4)*Xm(8)/(Xm(3)/2))/(Xm(1)*Xm(2)*Xm(9)*Xm(5))))+(1/((Zp(2)*Zp(4)*Zp(8)/(Zp(3)/2))/(Zp(1)*Zp(2)*Zp(9)*Zp(5)))));   %Xm to Zp
 Xm_Zm=1/((1/((Xm(2)*Xm(4)*Xm(8)/(Xm(3)/2))/(Xm(1)*Xm(2)*Xm(9)*Xm(5))))+(1/((Zm(2)*Zm(4)*Zm(8)/(Zm(3)/2))/(Zm(1)*Zm(2)*Zm(9)*Zm(5)))));   %Xm to Zm
 
 Zp_Yp=1/((1/((Zp(2)*Zp(3)*Zp(8)/(Zp(4)/2))/(Zp(1)*Zp(2)*Zp(9)*Zp(5))))+(1/((Yp(2)*Yp(4)*Yp(8)/(Yp(3)/2))/(Yp(1)*Yp(2)*Yp(9)*Yp(5)))));   %Zp to Yp
 Zm_Yp=1/((1/((Zm(2)*Zm(3)*Zm(8)/(Zm(4)/2))/(Zm(1)*Zm(2)*Zm(9)*Zm(5))))+(1/((Yp(2)*Yp(4)*Yp(8)/(Yp(3)/2))/(Yp(1)*Yp(2)*Yp(9)*Yp(5)))));   %Zm to Yp
 Zp_Ym=1/((1/((Zp(2)*Zp(3)*Zp(8)/(Zp(4)/2))/(Zp(1)*Zp(2)*Zp(9)*Zp(5))))+(1/((Ym(2)*Ym(4)*Ym(8)/(Ym(3)/2))/(Ym(1)*Ym(2)*Ym(9)*Ym(5)))));   %Zp to Ym
 Zm_Ym=1/((1/((Zm(2)*Zm(3)*Zm(8)/(Zm(4)/2))/(Zm(1)*Zm(2)*Zm(9)*Zm(5))))+(1/((Ym(2)*Ym(4)*Ym(8)/(Ym(3)/2))/(Ym(1)*Ym(2)*Ym(9)*Ym(5)))));   %Zm to Ym
 
 Zp_Sp=1/((1/((4*10^-4*Sp(8)/0.05)/((Sp(1)*Sp(2)*Sp(9)+0.45)*Sp(5))))+(1/((Sp(4)*Sp(2)*Sp(8)/(Sp(3)/2))/((Sp(1)*Sp(2)*Sp(9)+0.45)*Sp(5))))); %Zp to Sp
 Zm_Sp=1/((1/((4*10^-4*Sp(8)/0.05)/((Sp(1)*Sp(2)*Sp(9)+0.45)*Sp(5))))+(1/((Sp(4)*Sp(2)*Sp(8)/(Sp(3)/2))/((Sp(1)*Sp(2)*Sp(9)+0.45)*Sp(5))))); %Zm to Sp
 

 
 
 %% Worst case angles
for i=1:180

    
    
    for j=1:180
       A1(i,j)=Xp(1).*cosd(i).*cosd(j);
       A2(i,j)=Yp(1).*sind(i).*cosd(j);
       A3(i,j)=Zp(1).*cosd(i).*sind(j);
       Atot(i,j)=A1(i,j)+A2(i,j)+A3(i,j);
    end
    

end
[theta1, theta2] = find(ismember(Atot, max(Atot(:)))); % Theta1 is the angle between Xp and Y+, theta 2 is the angle between Xp and Zp 


%% Setup of vectors
 %temperature vector max
 Tmax_Xp=zeros(Period*n,1);
 Tmax_Xp(1)=T_i;
 %temperature vector min
 Tmin_Xp=zeros(Period*n,1);
 Tmin_Xp(1)=T_i;
 
  %temperature vector max
 Tmax_Xm=zeros(Period*n,1);
 Tmax_Xm(1)=T_i;
 %temperature vector min
 Tmin_Xm=zeros(Period*n,1);
 Tmin_Xm(1)=T_i;
 
  %temperature vector max
 Tmax_Yp=zeros(Period*n,1);
 Tmax_Yp(1)=T_i;
 %temperature vector min
 Tmin_Yp=zeros(Period*n,1);
 Tmin_Yp(1)=T_i;
 
  %temperature vector max
 Tmax_Ym=zeros(Period*n,1);
 Tmax_Ym(1)=T_i;
 %temperature vector min
 Tmin_Ym=zeros(Period*n,1);
 Tmin_Ym(1)=T_i;
 
  %temperature vector max
 Tmax_Zp=zeros(Period*n,1);
 Tmax_Zp(1)=T_i;
 %temperature vector min
 Tmin_Zp=zeros(Period*n,1);
 Tmin_Zp(1)=T_i;
 
  %temperature vector max
 Tmax_Zm=zeros(Period*n,1);
 Tmax_Zm(1)=T_i;
 %temperature vector min
 Tmin_Zm=zeros(Period*n,1);
 Tmin_Zm(1)=T_i;
 
 
   %temperature vector max
 Tmax_Sp=zeros(Period*n,1);
 Tmax_Sp(1)=T_i;
 %temperature vector min
 Tmin_Sp=zeros(Period*n,1);
 Tmin_Sp(1)=T_i;
 
 
 
 
 
 
 

 %% Incoming heat
%incoming heat day
 J_inc_max=(SF+EF+MF+R_M+R_E);
 J_inc_min=(SF+MF);
 %incoming heat night 
 Jn_inc_max=(EF+MF+R_E);
 Jn_inc_min=(MF);
 
 %max incoming radiation vector
 Inc_max=zeros(Period,1);
 Inc_max(1:Period-Eclipse_T)=J_inc_max;
 Inc_max(Period-Eclipse_T+1:Period)=Jn_inc_max;
 Inc_max=repmat(Inc_max,n,1)';
 
 %min incoming radiation vector
 Inc_min=zeros(Period,1);
 Inc_min(1:Period-Eclipse_T)=J_inc_min;
 Inc_min(Period-Eclipse_T+1:Period)=Jn_inc_min;
 Inc_min=repmat(Inc_min,n,1)';


%% Temperature simmulation

  k=1:(Period)*n-1;
 
 %Max temperature 
 for s=k
     Tmax_Xp(s+1)=Tmax_Xp(s)...
     +(((Inc_max(s)*Xp(1)*Xp(6)*cosd(theta1)*cosd(theta2))/(Xp(1)*Xp(2)*Xp(9)*Xp(5)))...
     -((sigma*Xp(7)*Xp(1)*Tmax_Xp(s)^4)/(Xp(1)*Xp(2)*Xp(9)*Xp(5)))...
     -((Tmax_Xp(s)-Tmax_Yp(s))*Xp_Yp)-((Tmax_Xp(s)-Tmax_Ym(s))*Xp_Ym)...
     -((Tmax_Xp(s)-Tmax_Zp(s))*Xp_Zp)-((Tmax_Xp(s)-Tmax_Zm(s))*Xp_Zm)...
     -((sigma*((Xp(7)+Xm(7))/2)*Xp(1)*(Tmax_Xp(s)^4-Tmax_Xm(s)^4))/(Xp(1)*Xp(2)*Xp(9)*Xp(5))));
 
     Tmax_Xm(s+1)=Tmax_Xm(s)...
     +(((0*Inc_max(s)*Xm(1)*Xm(6)*cosd(theta1)*cosd(theta2))/(Xm(1)*Xm(2)*Xm(9)*Xm(5)))...
     -((sigma*Xm(7)*Xm(1)*Tmax_Xm(s)^4)/(Xm(1)*Xm(2)*Xm(9)*Xm(5)))...
     -((Tmax_Xm(s)-Tmax_Yp(s))*Xm_Yp)-((Tmax_Xm(s)-Tmax_Ym(s))*Xm_Ym)...
     -((Tmax_Xm(s)-Tmax_Zp(s))*Xm_Zp)-((Tmax_Xm(s)-Tmax_Zm(s))*Xm_Zm)...
     -((sigma*((Xm(7)+Xp(7))/2)*Xm(1)*(Tmax_Xm(s)^4-Tmax_Xp(s)^4))/(Xm(1)*Xm(2)*Xm(9)*Xm(5))));
 
     Tmax_Ym(s+1)=Tmax_Ym(s)...
     +(((0*Inc_max(s)*Ym(1)*Ym(6)*cosd(theta1)*cosd(theta2))/(Ym(1)*Ym(2)*Ym(9)*Ym(5)))...
     -((sigma*Ym(7)*Ym(1)*Tmax_Ym(s)^4)/(Ym(1)*Ym(2)*Ym(9)*Ym(5)))...
     -((Tmax_Ym(s)-Tmax_Xp(s))*Xp_Ym)-((Tmax_Ym(s)-Tmax_Xm(s))*Xm_Ym)...
     -((Tmax_Ym(s)-Tmax_Zp(s))*Zp_Ym)-((Tmax_Ym(s)-Tmax_Zm(s))*Zm_Ym)...
     -((sigma*((Ym(7)+Yp(7))/2)*Ym(1)*(Tmax_Ym(s)^4-Tmax_Yp(s)^4))/(Ym(1)*Ym(2)*Ym(9)*Ym(5))));
 
     Tmax_Zm(s+1)=Tmax_Zm(s)...
     +(((0*Inc_max(s)*Zm(1)*Zm(6)*cosd(theta1)*cosd(theta2))/(Zm(1)*Zm(2)*Zm(9)*Zm(5)))...
     -((sigma*Zm(7)*Zm(1)*Tmax_Zm(s)^4)/(Zm(1)*Zm(2)*Zm(9)*Zm(5)))...
     -((Tmax_Zm(s)-Tmax_Xp(s))*Xp_Zm)-((Tmax_Zm(s)-Tmax_Xm(s))*Xm_Zm)...
     -((Tmax_Zm(s)-Tmax_Yp(s))*Zm_Yp)-((Tmax_Zm(s)-Tmax_Ym(s))*Zm_Ym)...
     -((sigma*((Zm(7)+Zp(7))/2)*Zm(1)*(Tmax_Zm(s)^4-Tmax_Zp(s)^4))/(Zm(1)*Zm(2)*Zm(9)*Zm(5))));
 
     Tmax_Zp(s+1)=Tmax_Zp(s)...
     +(((Inc_max(s)*Zp(1)*Zp(6)*cosd(theta1)*sind(theta2))/(Zp(1)*Zp(2)*Zp(9)*Zp(5)))...
     -((sigma*Zp(7)*Zp(1)*Tmax_Zp(s)^4)/(Zp(1)*Zp(2)*Zp(9)*Zp(5)))...
     -((Tmax_Zp(s)-Tmax_Xp(s))*Xp_Zp)-((Tmax_Zp(s)-Tmax_Xm(s))*Xm_Zp)...
     -((Tmax_Zp(s)-Tmax_Yp(s))*Zp_Yp)-((Tmax_Zp(s)-Tmax_Ym(s))*Zp_Ym)...
     -((sigma*((Zp(7)+Zm(7))/2)*Zp(1)*(Tmax_Zp(s)^4-Tmax_Zm(s)^4))/(Zp(1)*Zp(2)*Zp(9)*Zp(5))));
 
     Tmax_Yp(s+1)=Tmax_Yp(s)...
     +(((Inc_max(s)*Yp(1)*Yp(6)*sind(theta1)*cosd(theta2))/(Yp(1)*Yp(2)*Yp(9)*Yp(5)))...
     -((sigma*Yp(7)*Yp(1)*Tmax_Yp(s)^4)/(Yp(1)*Yp(2)*Yp(9)*Yp(5)))...
     -((Tmax_Yp(s)-Tmax_Xp(s))*Xp_Yp)-((Tmax_Yp(s)-Tmax_Xm(s))*Xm_Yp)...
     -((Tmax_Yp(s)-Tmax_Zp(s))*Zp_Yp)-((Tmax_Yp(s)-Tmax_Zm(s))*Zm_Yp)...
     -((sigma*((Yp(7)+Ym(7))/2)*Yp(1)*(Tmax_Yp(s)^4-Tmax_Ym(s)^4))/(Yp(1)*Yp(2)*Yp(9)*Yp(5))));
 
     Tmax_Sp(s+1)=Tmax_Sp(s)...
     +(((Inc_max(s)*Sp(1)*Sp(6)*cosd(panel_incl))/((Sp(1)*Sp(2)*Sp(9)+0.45)*Sp(5)))...
     -(((sigma*Sp(7)*Sp(1)*Tmax_Sp(s)^4)+(sigma*Sp(7)*Sp(1)*Tmax_Sp(s)^4))/((Sp(1)*Sp(2)*Sp(9)+0.45)*Sp(5)))...
     -((Tmax_Sp(s)-Tmax_Zp(s))*Zp_Sp)...
     -((0*sigma*((Yp(7)+Ym(7))/2)*Yp(1)*(Tmax_Yp(s)^4-Tmax_Ym(s)^4))/(Yp(1)*Yp(2)*Yp(9)*Yp(5))));

 
 
 
 
 
 
 end
    
 
 %min temperature
 
%   for s=k
%      T_min(s+1)=T_min(s)+((Inc_min(s)-sigma.*(T_min(s)^4).*(T_emi)+IH)./(Mass.*SPH));
%   end


%% Plotting 

time=1:Period*n;
%  figure
 
%  title('Temperature plots over 10 orbits')
%  x=plot(time,Tmax_Xp,time,Tmax_Xm,time,Tmax_Yp,time,Tmax_Ym,time,Tmax_Zp,time,Tmax_Zm);
%  title('Temperature plots over 10 orbits')
%  xlabel('Time after initial state [s]')
%  ylabel('Temperature [K]')
%  x(1).LineWidth = 1;
%  x(4).LineWidth=1;
%  x(5).LineWidth=1;
%  x(2).LineWidth=1;
%  x(3).LineWidth=1;
%  grid on
%  legend('Xp panel','Xm panel','Yp panel','Ym panel','Zp panel','Zm panel')

%  figure
%  y=plot(time,Tmax_Sp);
%  title('Temperature plots over 10 orbits')
%  xlabel('Time after initial state [s]')
%  ylabel('Temperature [K]')
% 
%  grid on
%  legend('Sp panel','Sm panel')
 figure;
 Temps_out = Tmax_Sp(Period+1:Period*2);
 
 
plot(1:Period,Temps_out(1:end))
 
 

end

