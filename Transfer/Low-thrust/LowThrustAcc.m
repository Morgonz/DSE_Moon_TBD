function [f] = LowThrustAcc(t,state)

xp = state(1); %x-position of the satellite
yp = state(2); %y-position of the satellite
zp = state(3); %z-position of the satellite
xv = state(4); %speed in x-direction of the satellite
yv = state(5); %speed in y-direction of the satellite
zv = state(6); %speed in z-direction of the satellite

xM = state(7); %x-position of the Moon
yM = state(8); %y-position of the Moon
zM = state(9); %z-position of the Moon
xvM = state(10); %speed in x-direction of the Moon
yvM = state(11); %speed in y-direction of the Moon
zvM = state(12); %speed in z-direction of the Moon

Mu_Earth = 3.98574405E+14; %m^3 s^-2 Gravitational parameter of the Earth
Mu_Moon = 4.902801e12; %m^3 s^-2 Gravitational parameter of the Moon

v = sqrt(xv^2 + yv^2 + zv^2); %velocity satellite
pS = sqrt(xp^2 + yp^2 + zp^2); %position satellite
pM = sqrt(xM^2 + yM^2 + zM^2); %position of the Moon
dSM = sqrt((xM-xp)^2 + (yM-yp)^2 + (zM-zp)^2); %relative distance satellite to the moon

asE = -(Mu_Earth/(pS*pS*pS)); %acceleration of the satellite due to the Earth
asM = -(Mu_Moon/(dSM^3)); %acceleration of the satellite due to the Moon
aM  = -((Mu_Earth + Mu_Moon)/(pM^3)); %acceleration of the Moon

%mass calculations
Mdry = 2000; %kg dry mass estimation
Mprop = 1000; %kg propellant mass estimation
m = 1000; %kg initial mass

%thrust force
N = 7*0.09; %N PPS-1350 Hall effect 7 thrusters for same dry mass-thruster ratio as SMART-1
aN = N/m; %acceleration caused by the thrust force

xa = asE.* xp + asM.*(xp-xM) + aN*(xv/v); %accelaration in x-direction of the satellite
ya = asE.* yp + asM.*(yp-yM) + aN*(yv/v); %acceleration in y-direction of the satellite
za = asE.* zp + asM.*(zp-zM) + aN*(zv/v); %acceleration in z-direction of the satellite

xaM = aM*xM; %acceleration in x-direction of the Moon
yaM = aM*yM; %acceleration in y-direction of the Moon
zaM = aM*zM; %acceleration in z-direction of the Moon


f = [xv yv zv xa ya za xvM yvM zvM xaM yaM zaM]';