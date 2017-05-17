tic
time = [];
for theta = linspace(-180,180,3601)
    Mu_Earth = 3.98574405E+14; %m^3 s^-2
    Mu_Moon = 4.902801e12; %m^3 s^-2
    R_Earth = 6371000; %m
    R_Moon = 1737000; %m

%Circular Orbit initial conditions
    h = 500000; %m
    V0 = sqrt(Mu_Earth/(R_Earth + h)); %m/s
    r1 = R_Earth + h; %m

    %orbit of the Moon
    rM = 385000600; %m
    vM = sqrt(Mu_Earth/rM); %m/s speed of rotation of the Moon

    %transfer orbit
    aT  = (r1+rM)/2; 
    DV1 = sqrt(Mu_Earth/r1)*(sqrt((2*rM)/(r1+rM))-1); %First Delta V
    T = 864000; %s time of simulation
    V1 = V0+DV1; %velocity just after the initial burn

    %orbit around the Moon
    h1 = 1000000; %m
    r3 = h1 + R_Moon; %radius of the orbit around the Moon
    V3 = sqrt(Mu_Moon/r3); %m/s orbiting speed in Moon reference frame


    %options of the integrator
    options1 = odeset('RelTol', 1e-12, 'Events', @CrossMoonOrbit);


    [t2,y2,te] = ode45(@HohmannAcc,[0 T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V1*sin(theta*(pi/180)) V1*cos(theta*(pi/180)) 0 rM 0 0 0 vM 0],options1); %Transfer Orbit
    if te>1000
        xp2 = y2(length(y2),1); %position in x-direction as transfer orbit crosses target orbit
        yp2 = y2(length(y2),2); %position in y-direction as transfer orbit crosses target orbit
        zp2 = y2(length(y2),3); %position in z-direction as transfer orbit crosses target orbit
        xV2 = y2(length(y2),4); %velocity in x direction as transfer orbit crosses target orbit
        yV2 = y2(length(y2),5); %velocity in y direction as transfer orbit crosses target orbit
        zV2 = y2(length(y2),6); %velocity in z direction as transfer orbit crosses target orbit
        xM2 = y2(length(y2),7); %x-coordinate of Moon
        yM2 = y2(length(y2),8); %y-coordinate of Moon
        zM2 = y2(length(y2),9); %z-coordinate of Moon
        xVM2 = y2(length(y2),10); %velocity in x-direction of the Moon
        yVM2 = y2(length(y2),11); %velocity in y-direction of the Moon
        zVM2 = y2(length(y2),12); %velocity in z-direction of the Moon
        
        if xp2>xM2
            theta2=atan((yp2-yM2)/(xp2-xM2));
        else
            theta2=atan((yp2-yM2)/(xp2-xM2))+pi;
        end
        xV3 = V3*sin(theta2); %target velocity in x-direction in Moon reference frame
        yV3 = -V3*cos(theta2); %target velocity in y-direction in Moon reference frame
        xVtar = xV3+xVM2; %target velocity in x-direction in Earth reference frame
        yVtar = yV3+yVM2; %target velocity in y-direction in Earth reference frame
        zVtar = 0;
        
        xDV2 = xVtar - xV2; %delta V in x-direction
        yDV2 = yVtar - yV2; %delta V in y-direction
        zDV2 = zVtar - zV2; %delta V in z-direction
        
        DV2 = sqrt(xDV2^2 + yDV2^2 + zDV2^2);
        
        time = [time; [theta DV2 xDV2 yDV2 zDV2]];
        
         
    end
end
toc
    %Define Event Function, target orbit at 1000 km above Moon surface
    function [value,isterminal,direction] = CrossMoonOrbit(t2,y2)
    value = sqrt((y2(1)-y2(7))^2 + (y2(2)-y2(8))^2 + (y2(3)-y2(9))^2)-2737000;
    isterminal = 1;
    direction = 1;
    end