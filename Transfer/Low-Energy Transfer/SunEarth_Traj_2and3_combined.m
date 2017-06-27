clear all;
close all;
clc;

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Sun = 1.327124E+20; %m^3 s^-2 from wikipedia
R_Earth = 6371000; %m
R_Sun = 695508000; %m https://solarsystem.nasa.gov/planets/sun/facts

%Earth initial conditions
rE = 149.60E+9; %m radius of the Earth https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
vE = sqrt((Mu_Earth+Mu_Sun)/rE); %rotational speed of the Earth

%Sun-Earth L2
rL = rE*(Mu_Earth/(3*Mu_Sun))^(1/3);

% SUN
[X,Y,Z] = sphere(40);

%options of the integrator
options1 = odeset('RelTol', 2.22045e-14);
options2 = odeset('RelTol', 2.22045e-14, 'Events', @toofar);

%define delta v of kick adn total time
dV_kick = -10 %in delta V
%T_days = 183;
T_days = 181.1918
T = 60*60*24*T_days; %s


%for t_manouvre = linspace(1,180,180)
for t_manouvre = [ 59.14 ]
    
    t = 60*60*24*t_manouvre; %s

    [t1,y1] = ode113(@SunEarthAcc, [0 t], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options1);
    yrot1 = RotatingFrameSunEarth(y1);
    
    
    
    % calculating velocities, use only the X and Y velocities
    V_abs = sqrt(y1(end,10)^2+y1(end,11)^2+y1(end,11)^2);
    V_abs_xy = sqrt(y1(end,10)^2+y1(end,11)^2);
    dv_x = (y1(end,10)/V_abs_xy)*dV_kick;
    dv_y = (y1(end,11)/V_abs_xy)*dV_kick;
    dv_z = (y1(end,12)/V_abs_xy)*dV_kick;
    
    %without z-velocity 
    [t2,y2] = ode113(@SunEarthAcc, [0 T-t], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+dv_x y1(end,11)+dv_y y1(end,12)],options1);
    %with z-velocity included
    %[t2,y2] = ode113(@SunEarthAcc, [0 T-t], [y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5) y1(end,6) y1(end,7) y1(end,8) y1(end,9) y1(end,10)+dv_x y1(end,11)+dv_y y1(end,12)+dv_z],options2);
    yrot2 = RotatingFrameSunEarth(y2);
    
    if yrot2(end,7)<1.514E+11 & yrot2(end,8)> 2.0E+8 & yrot2(end,8)< 5.0E+8
        hold on
        Sat_trajectory = plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9),'r-','DisplayName','Satellite trajectory');
        burn4 = plot3(yrot2(end,7),yrot2(end,8),yrot2(end,9),'r<','MarkerFaceColor','r','DisplayName','Matching burn with trajectory part 3');
        burn3 = plot3(yrot1(end,7),yrot1(end,8),yrot1(end,9),'k<','MarkerFaceColor','k','DisplayName','Burn to get out of the halo L2 orbit');
        %plot3(yrot1(end,7),yrot1(end,8),yrot1(end,9),'k<','DisplayName','orbit kick');
        %plot3(yrot2(:,7),yrot2(:,8),yrot2(:,9),'r-','DisplayName',['manouvre at', num2str(t_manouvre) ]);
    end
end
%figure
%plot(t,yrot2(:,4))
%figure
%plot(t,yrot2(:,5))

%only the Halo orbit itself
T_halo = 178*60*60*24
[t_halo,y_halo] = ode113(@SunEarthAcc, [0 T_halo], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options1);
yrot_halo = RotatingFrameSunEarth(y_halo);
halo_orbit = plot3(yrot_halo(:,7),yrot_halo(:,8),yrot_halo(:,9),'b--','DisplayName','Sun-Earth L2 halo orbit');

%% For total plot
%Trajectory part 1
y_inertial_bart = importdata('State_Earth_to_L2_PILSBAAS_LAURENCE.mat')
yrot_bart = RotatingFrameSunEarth(y_inertial_bart);
%plot3(yrot_bart(:,7),yrot_bart(:,8),yrot_bart(:,9),'r-','DisplayName','Trajectory part 1');
plot3(yrot_bart(:,7),yrot_bart(:,8),yrot_bart(:,9),'r-');
burn2 = plot3(yrot_bart(1,7),yrot_bart(1,8),yrot_bart(1,9),'<','Color',[0.0 0.8 1],'MarkerFaceColor',[0.0 0.8 1], 'DisplayName','Insertion burn into the halo L2 orbit');
burn1 = plot3(yrot_bart(end,7),yrot_bart(end,8),yrot_bart(end,9),'g<','MarkerFaceColor','g','DisplayName','LEO to Sun-Earth L2 halo burn');

%Trajecoroy intermediate part in halo L2
T_intermediate_begin = 5*60*60*24
T_intermediate_end = 59.14*60*60*24
T_intermediate = T_intermediate_end-T_intermediate_begin
[t_int1,y_int1] = ode113(@SunEarthAcc, [0 T_intermediate_begin], [rE 0 0 0 vE 0 rE+rL-120000000 0 0 1.06749 vE+426 1.260005700065],options1);
[t_int2,y2_int2] = ode113(@SunEarthAcc, [0 T_intermediate], [y_int1(end,1) y_int1(end,2) y_int1(end,3) y_int1(end,4) y_int1(end,5) y_int1(end,6) y_int1(end,7) y_int1(end,8) y_int1(end,9) y_int1(end,10) y_int1(end,11) y_int1(end,12)],options1);
yrot_int = RotatingFrameSunEarth(y2_int2);
plot3(yrot_int(:,7),yrot_int(:,8),yrot_int(:,9),'r-','DisplayName','Intermediate part in halo');
%burn2 = plot3(yrot_int(1,7),yrot_int(10,8),yrot_int(10,9),'<','Color',[0.0 0.8 1],'MarkerFaceColor',[0.0 0.8 1], 'DisplayName','Insertion burn into the halo L2 ');

%%
%Trajectory part 3
load('ZIEKE_FIGURES_VOOR_PILSBAAS_LAURENCE_DIE_EEN_HOOFD_HEEFT_DIE_NIET_ALTIJD_MEE_WERKT.mat')
plot3(y2_p3_new(:,7),y2_p3_new(:,8),y2_p3_new(:,9),'r-')
burn5 = plot3(y2_p3_new(1,7),y2_p3_new(1,8),y2_p3_new(1,9),'<','Color',[0.8 0. 1],'MarkerFaceColor',[0.8 0.0 1], 'DisplayName','Insertion burn into the Earth-Moon halo L1 orbit')
plot3(y2_p3_new(:,1),y2_p3_new(:,2),y2_p3_new(:,3),'k-')
%plot3()




%%

%title(['Asymptotic trajectories leaving Earth-Sun L2 with an deltaV kick of: ' num2str(dV_kick) 'm/s'])
title('Trajectory part 1 and 2 combined in Sun-centered rotating frame, including burns')



Earth_yaxis = plot([1.496e+11 1.496e+11],[0 5e+8],'k--','DisplayName',' Earth inertial y-axis')


Earth = plot3(rE,0,0,'b.','MarkerSize',25,'DisplayName','Earth location')
L2 = plot3(rE+rL,0,0,'k*','DisplayName','Sun-Earth L2')
%surf(X*R_Sun,Y*R_Sun, Z*R_Sun,'DisplayName','Sun position')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
%legend([L2 Earth],{'Earth-Sun L2','Earth'})
legend([Sat_trajectory halo_orbit Earth_yaxis L2 Earth burn1 burn2 burn3 burn4 burn5])
%legend('show')
axis equal
axis vis3d



function [value,isterminal,direction] = toofar(t,y)
%disp('value found is')
if y(1)>0
    tet = atan(y(2)/y(1));
elseif y(1)<0
    tet = atan(y(2)/y(1))+pi;
else
    if y(2)>0
        tet = pi/2;
    else
        tet = -pi/2;   
    end
end
xEr1 = y(1)*cos(-tet)-y(2)*sin(-tet);
yEr1 = y(1)*sin(-tet)+y(2)*cos(-tet);
xSr1 = y(7)*cos(-tet)-y(8)*sin(-tet);
ySr1 = y(7)*sin(-tet)+y(8)*cos(-tet);

%value = 2*1.496547398746715e+09 - sqrt((xSr1-1.5110E+11)^2 + (ySr1)^2);
value = 1.496547398746715e+09 - abs(xSr1-1.5110E+11);
%value =  xSr1-1.49600e+11

isterminal = 1;
direction = 0;
end


