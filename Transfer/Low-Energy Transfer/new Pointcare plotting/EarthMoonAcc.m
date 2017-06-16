function [f] = EarthMoon(t, state)
Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Moon = 4.902801e12; %m^3 s^-2

xM = state(1);
yM = state(2);
zM = state(3);
xvM = state(4);
yvM = state(5);
zvM = state(6);
xS = state(7);
yS = state(8);
zS = state(9);
xvS = state(10);
yvS = state(11);
zvS = state(12);

vM = sqrt(xvM^2 + yvM^2 + zvM^2); %velocity Moon, not necesserary to calculate actually
pM = sqrt(xM^2 + yM^2 + zM^2); %Distance Moon-Earth
vS = sqrt(xvS^2 + yvS^2 + zvS^2); %velocity of the satellite
pS = sqrt(xS^2 + yS^2 + zS^2); %Distance Satellite-Earth
pSM = sqrt((xS-xM)^2 + (yS-yM)^2 + (zS-zM)^2); %distance Satellite-Moon

aME = -((Mu_Earth+Mu_Moon)/(pM^3)); %acceleration of the Moon due to the Earth
aSE = -((Mu_Earth)/(pS^3)); %acceleration of the satellite due to the Earth
aSM = -((Mu_Moon)/(pSM^3)); %acceleration of the satellite due to the Moon

xaM = aME.* xM ; %accelaration in x-direction of the Moon
yaM = aME.* yM ; %acceleration in y-direction of the Moon
zaM = aME.* zM ; %acceleration in z-direction of the Moon
xaS = aSE.* xS + aSM*(xS-xM); %acceleration in x-direction of the satellite
yaS = aSE.* yS + aSM*(yS-yM); %acceleration in y-direction of the satellite
zaS = aSE.* zS + aSM*(zS-zM); %acceleration in z-direction of the satellite

f = [xvM yvM zvM xaM yaM zaM xvS yvS zvS xaS yaS zaS]';
