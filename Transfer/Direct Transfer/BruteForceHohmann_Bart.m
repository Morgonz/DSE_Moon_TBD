clear all;
close all;
clc;

time = [];
for DV1 = 3059.93304058853%linspace(3059.93304058853-7,3059.93304058853+12,20)
    time1 = [];
    for theta = -122.6%linspace(-180,180,3601)

        Mu_Earth = 3.98574405E+14; %m^3 s^-2
        Mu_Moon = 4.902801e12; %m^3 s^-2
        R_Earth = 6371000; %m
        R_Moon = 1737000; %m

    %Circular Orbit initial conditions
        h = 500000; %m
        V0 = sqrt(Mu_Earth/(R_Earth + h)); %m/s
        r1 = R_Earth + h; %m

        %orbit of the Moon
        %orbit of the Moon
        aM = 384400000; %m
        eM = 0.0554; %Eccentricity
        rMa = aM*(1+eM);%Apocentre
        rMp = aM*(1-eM); %Pericentre
        vM = sqrt(((2*(Mu_Earth+Mu_Moon))/rMa)-((Mu_Earth+Mu_Moon)/aM)); %m/s speed of rotation of the Moon at the pericentre
        
        %transfer orbit
        %aT  = (r1+rM)/2; 
        %DV1 = sqrt(Mu_Earth/r1)*(sqrt((2*rM)/(r1+rM))-1); %First Delta V
        T = 864000; %s time of simulation
        V1 = V0+DV1; %velocity just after the initial burn

        %Frozen orbit around the Moon
        h3 = 1629000; %km frozen orbit
        r3 = h3+R_Moon;
        V3 = sqrt(Mu_Moon/r3);
        i_change = (46)*(pi/180);


        %options of the integrator
        options1 = odeset('RelTol', 1e-12, 'Events', @CrossMoonOrbit);


        [t2,y2,te] = ode45(@HohmannAcc,[0 T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V1*sin(theta*(pi/180)) V1*cos(theta*(pi/180)) 0 rMa 0 0 0 vM 0],options1); %Transfer Orbit
        if te>1000
            xp2 = y2(end,1); %position in x-direction as transfer orbit crosses target orbit
            yp2 = y2(end,2); %position in y-direction as transfer orbit crosses target orbit
            zp2 = y2(end,3); %position in z-direction as transfer orbit crosses target orbit
            xV2 = y2(end,4); %velocity in x direction as transfer orbit crosses target orbit
            yV2 = y2(end,5); %velocity in y direction as transfer orbit crosses target orbit
            zV2 = y2(end,6); %velocity in z direction as transfer orbit crosses target orbit
            xM2 = y2(end,7); %x-coordinate of Moon
            yM2 = y2(end,8); %y-coordinate of Moon
            zM2 = y2(end,9); %z-coordinate of Moon
            xVM2 = y2(end,10); %velocity in x-direction of the Moon
            yVM2 = y2(end,11); %velocity in y-direction of the Moon
            zVM2 = y2(end,12); %velocity in z-direction of the Moon

            if xp2>xM2
                theta2=atan((yp2-yM2)/(xp2-xM2));
            else
                theta2=atan((yp2-yM2)/(xp2-xM2))+pi;
            end
            zV3 =V3*sin(i_change);%target velocity in z-direction in Moon rference frame
            V3b =V3*cos(i_change);%target velocity in x-y plane
            xV3 = -V3b*sin(theta2); %target velocity in x-direction in Moon reference frame
            yV3 = V3b*cos(theta2); %target velocity in y-direction in Moon reference frame

            xVtar = xV3+xVM2; %target velocity in x-direction in Earth reference frame
            yVtar = yV3+yVM2; %target velocity in y-direction in Earth reference frame
            zVtar = zV3;

            xDV2 = xVtar - xV2; %delta V in x-direction
            yDV2 = yVtar - yV2; %delta V in y-direction
            zDV2 = zVtar - zV2; %delta V in z-direction

            DV2 = sqrt(xDV2^2 + yDV2^2 + zDV2^2);

            time = [time; [DV1 theta DV2 xDV2 yDV2 zDV2]];
            time1 = [time1; [DV1 theta DV2]];
        else
            time1 = [time1; [DV1 theta NaN]];
        end
        disp(theta), disp(DV1)
    end
    hold on
    plot(time1(:,2),time1(:,3),'DisplayName',num2str(round(time1(1,1),1)))
    lekkur = legend('show');
    lekkurt = get(lekkur,'Title');
    set(lekkurt,'String','\Delta V_{21}')
end


    %Define Event Function, target orbit above Moon surface
    function [value,isterminal,direction] = CrossMoonOrbit(t2,y2)
    value = sqrt((y2(1)-y2(7))^2 + (y2(2)-y2(8))^2 + (y2(3)-y2(9))^2)-3366000;
    isterminal = 1;
    direction = 0;
    end
