clear all;
close all;
clc;

Mu_Moon = 4.902801e12; %m^3 s^-2
R_Moon = 1737000; %m

%initial orbit
h = 1629000; %m
r1 = h+R_Moon; %m
v1 = sqrt(Mu_Moon/r1); %m/s
a = 3271838.37;
vphase = sqrt((2*Mu_Moon/r1)-(Mu_Moon/a));

T = 167936.6544; %s
Tcirc = 2*pi*sqrt((r1^3)/Mu_Moon);
options = odeset('RelTol', 1e-12);

[t,y] = ode113(@PhaseAcc, [0 T], [r1 0 0 0 v1 0],options);
[t2,y2] = ode113(@PhaseAcc, [0 T], [r1 0 0 0 v1-(v1-vphase) 0],options);

phase1 = [];
for i = linspace(1,length(y),length(y))
    if y(i,1)>0 && y(i,2)>0
        phase = (180/pi)*atan(y(i,2)/y(i,1)); %in degrees
    elseif y(i,1)==0 && y(i,2)>0
        phase = 90;
    elseif y(i,1)<0 && y(i,2)>0
        phase = 180+((180/pi)*atan(y(i,2)/y(i,1))); %in degrees
    elseif y(i,2)==0 && y(i,1)<0
        phase = 180;
    elseif y(i,1)<0 && y(i,2)<0
        phase = 180+((180/pi)*atan(y(i,2)/y(i,1))); %in degrees
    elseif y(i,1)==0 && y(i,2)<0
        phase = 270;
    elseif y(i,1)>0 && y(i,2)<0
        phase = 360+((180/pi)*atan(y(i,2)/y(i,1)));
    elseif y(i,1)>0 && y(i,2)==0
        phase = 0;
    end
    phase1 = [phase1; phase];
end
phase2 = [];
for i = linspace(1,length(y2),length(y2))
    if y2(i,1)>0 && y2(i,2)>0
        phase = (180/pi)*atan(y2(i,2)/y2(i,1)); %in degrees
    elseif y2(i,1)==0 && y2(i,2)>0
        phase = 90;
    elseif y2(i,1)<0 && y2(i,2)>0
        phase = 180+((180/pi)*atan(y2(i,2)/y2(i,1))); %in degrees
    elseif y2(i,2)==0 && y2(i,1)<0
        phase = 180;
    elseif y2(i,1)<0 && y2(i,2)<0
        phase = 180+((180/pi)*atan(y2(i,2)/y2(i,1))); %in degrees
    elseif y2(i,1)==0 && y2(i,2)<0
        phase = 270;
    elseif y2(i,1)>0 && y2(i,2)<0
        phase = 360+((180/pi)*atan(y2(i,2)/y2(i,1)));
    elseif y2(i,1)>0 && y2(i,2)==0
        phase = 360;
    end
    phase2 = [phase2; phase];
end

tphase = 10*t/T;
tphase2 = 10*t2/T;

pcirc = (360/Tcirc)*t2;
for i=linspace(1,length(pcirc),length(pcirc))
    while pcirc(i)>360
        pcirc(i)=pcirc(i)-360;
    end
end

pdif = phase2-pcirc;
for i=linspace(1,length(pdif),length(pdif))
    if abs(pdif(i))>2 && pdif(i)<0
        pdif(i)=pdif(i)+360;
    end
end
[X,Y,Z] = sphere(40);

hold on
b1=plot(y(:,1),y(:,2),'c');
b2=plot(y2(:,1),y2(:,2),'b');
b3=plot(y(end,1),y(end,2),'oc');
b4=plot(y2(end,1),y2(end,2),'*b');
b5=plot(r1,0,'<k');
b6=surf(X*R_Moon,Y*R_Moon,Z*R_Moon);
axis equal
axis vis3d
hold off
xlabel('x [m]')
ylabel('y [m]')
legend([b1,b2,b3,b4,b5],{'Circular orbit', 'Phasing orbit','End location circular orbit','End location phasing orbit','Burns'})
title('Phasing manoeuvre')

% figure
% hold on
% b1=plot(tphase,phase1,'c');
% b2=plot(tphase2,phase2);
% b3=plot(tphase2,pdif);
% xlim([0.1,tphase2(end)])
% hline = refline([0 150]);
% set(hline,'LineStyle','--')
% ylim([0 450]);
% xlabel('Number of periods in phasing orbit [-]')
% ylabel('Phase [\circ]')
% title('Phasing manoeuvre')
% legend([b1,b2,b3,hline],{'Phase in circular orbit', 'Phase in phasing orbit', 'Phase difference', 'Desired phase difference'});

