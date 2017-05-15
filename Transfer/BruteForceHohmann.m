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
    v3 = sqrt(Mu_Moon/r3); %m/s orbiting speed in Moon reference frame


    %options of the integrator
    options1 = odeset('RelTol', 1e-12, 'Events', @CrossMoonOrbit);


    [t2,y2,te] = ode45(@HohmannAcc,[0 T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V1*sin(theta*(pi/180)) V1*cos(theta*(pi/180)) 0 rM 0 0 0 vM 0],options1); %Transfer Orbit
    if te>1000
        x2 = y2(length(y2),4);
        y2 = y2(length(y2),4);
        xV2 = y2(length(y2),4);
        xV2 = y2(length(y2),4); %velocity in x direction as transfer orbit crosses target orbit
        yV2 = y2(length(y2),5); %velocity in y direction as transfer orbit crosses target orbit
        zV2 = y2(length(y2),6); %velocity in z direction as transfer orbit crosses target orbit
        
    end
end

    %Define Event Function, target orbit at 1000 km above Moon surface
    function [value,isterminal,direction] = CrossMoonOrbit(t2,y2)
    value = sqrt((y2(1)-y2(7))^2 + (y2(2)-y2(8))^2 + (y2(3)-y2(9))^2)-2737000;
    isterminal = 1;
    direction = 1;
    end