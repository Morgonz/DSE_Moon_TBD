function [f] = SunEarth(t, state)
Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Sun = 1.327124E+20; %m^3 s^-2 from wikipedia
R_Earth = 6371000; %m
R_Sun = 695508000; %m https://solarsystem.nasa.gov/planets/sun/facts

xE = state(1);
yE = state(2);
zE = state(3);
xvE = state(4);
yvE = state(5);
zvE = state(6);
xS = state(7);
yS = state(8);
zS = state(9);
xvS = state(10);
yvS = state(11);
zvS = state(12);

vE = sqrt(xvE^2 + yvE^2 + zvE^2); %velocity Earth
pE = sqrt(xE^2 + yE^2 + zE^2); %Distance Earth-Sun
vS = sqrt(xvS^2 + yvS^2 + zvS^2); %velocity of the satellite
pS = sqrt(xS^2 + yS^2 + zS^2); %Distance Satellite-Sun
pSE = sqrt((xS-xE)^2 + (yS-yE)^2 + (zS-zE)^2); %distance Satellite-Earth

aE = -((Mu_Sun+Mu_Earth)/(pE*pE*pE)); %acceleration of the Earth due to the Sun
aSS = -((Mu_Sun)/(pS*pS*pS)); %acceleration of the satellite due to the Sun
aSE = -((Mu_Earth)/(pSE*pSE*pSE)); %acceleration of the satellite due to the Earth

xaE = aE.* xE ; %accelaration in x-direction of the Earth
yaE = aE.* yE ; %acceleration in y-direction of the Earth
zaE = aE.* zE ; %acceleration in z-direction of the Earth
xaS = aSS.* xS + aSE*(xS-xE); %acceleration in  x-direction of the satellite
yaS = aSS.* yS + aSE*(yS-yE); %acceleration in y-direction of the satellite
zaS = aSS.* zS + aSE*(zS-zE); %acceleration in z-direction of the satellite

f = [xvE yvE zvE xaE yaE zaE xvS yvS zvS xaS yaS zaS]';
