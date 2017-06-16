function [f] = PhaseAcc(t,state)

xp = state(1);
yp = state(2);
zp = state(3);
vx = state(4);
vy = state(5);
vz = state(6);

Mu_Moon = 4.902801e12; %m^3 s^-2

p = sqrt(xp^2 + yp^2 + zp^2);
asM = -(Mu_Moon/(p^3)); %acceleration of the satellite due to the Moon

xa = asM* xp; %accelaration in x-direction of the satellite
ya = asM* yp; %acceleration in y-direction of the satellite
za = asM* zp; %acceleration in z-direction of the satellite
f = [vx vy vz xa ya za]';
