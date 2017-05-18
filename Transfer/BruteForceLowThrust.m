time = [];
Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Moon = 4.902801e12; %m^3 s^-2
R_Earth = 6371000; %m
R_Moon = 1737000; %m
for theta = linspace(-180,180,37)
    %Circular Orbit initial conditions
    h = 500000; %m
    V0 = sqrt(Mu_Earth/(R_Earth + h)); %m/s
    r1 = R_Earth + h; %m
    
    %orbit of the Moon
    rM = 385000600; %m
    vM = sqrt(Mu_Earth/rM); %m/s speed of rotation of the Moon

    %options of the integrator
    options2 = odeset('RelTol', 1e-9, 'Events', @CrossMoonOrbit);
    T = 10*8640000; %s time of simulation

    [t2,y2,te] = ode45(@LowThrustAcc,[0,T],[r1*cos(theta*(pi/180)) r1*sin(theta*(pi/180)) 0 -V0*sin(theta*(pi/180)) V0*cos(theta*(pi/180))  0 rM 0 0 0 vM 0],options2); %Transfer orbit

    if te>1000
        time = [time; theta];
    end
    disp(theta)
end
%Define Event Function
function [value,isterminal,direction] = CrossMoonOrbit(t2,y2)
value = sqrt((y2(1)-y2(7))^2 + (y2(2)-y2(8))^2 + (y2(3)-y2(9))^2)-1071000; %stop burning when entering sphere of influence of the Moon
isterminal = 1;
direction = 1;
end
