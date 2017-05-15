function [f] = HohmannAcc(t,state)

xp = state(1);
yp = state(2);
zp = state(3);
xv = state(4);
yv = state(5);
zv = state(6);

xM = state(7);
yM = state(8);
zM = state(9);
xvM = state(10);
yvM = state(11);
zvM = state(12);

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Moon = 4.902801e12; %m^3 s^-2

v = sqrt(xv^2 + yv^2 + zv^2); %velocity satellite
pS = sqrt(xp^2 + yp^2 + zp^2); %position satellite
pM = sqrt(xM^2 + yM^2 + zM^2); %position of the Moon
dSM = sqrt((xM-xp)^2 + (yM-yp)^2 + (zM-zp)^2); %relative distance satellite to the moon

asE = -(Mu_Earth/(pS*pS*pS)); %acceleration of the satellite due to the Earth
asM = -(Mu_Moon/(dSM^3)); %acceleration of the satellite due to the Moon
aM  = -((Mu_Earth + Mu_Moon)/(pM^3)); %acceleration of the Moon

xa = asE.* xp + asM.*(xp-xM); %accelaration in x-direction of the satellite
ya = asE.* yp + asM.*(yp-yM); %acceleration in y-direction of the satellite
za = asE.* zp + asM.*(zp-zM); %acceleration in z-direction of the satellite

xaM = aM*xM; %acceleration in x-direction of the Moon
yaM = aM*yM; %acceleration in y-direction of the Moon
zaM = aM*zM; %acceleration in z-direction of the Moon

f = [xv yv zv xa ya za xvM yvM zvM xaM yaM zaM]';