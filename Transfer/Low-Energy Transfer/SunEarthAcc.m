function [f] = SunEarthAcc(t,state)

xE = state(1);
yE = state(2);
zE = state(3);
xvE = state(4);
yvE = state(5);
zvE = state(6);
xM = state(7);
yM = state(8);
zM = state(9);
xvM = state(10);
yvM = state(11);
zvM = state(12);
xS = state(13);
yS = state(14);
zS = state(15);
xvS = state(16);
yvS = state(17);
zvS = state(18);

Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Sun = 1.327124E+20; %m^3 s^-2 from wikipedia
Mu_Moon = 4.902801e12; %m^3 s^-2


vE = sqrt(xvE^2 + yvE^2 + zvE^2); %velocity Earth
pE = sqrt(xE^2 + yE^2 + zE^2); %Distance Earth-Sun
vM = sqrt(xvM^2 + yvM^2 +zvM^2); %velocity of the Moon
pM = sqrt(xM^2 + yM^2 +zM^2); %Distance Moon-Sun
vS = sqrt(xvS^2 + yvS^2 + zvS^2); %velocity of the satellite
pS = sqrt(xS^2 + yS^2 + zS^2); %Distance Satellite-Sun
pSE = sqrt((xS-xE)^2 + (yS-yE)^2 + (zS-zE)^2); %distance Satellite-Earth
pSM = sqrt((xS-xM)^2 + (yS-yM)^2 + (zS-zM)^2); %distance Satellite-Moon
pEM = sqrt((xE-xM)^2 + (yE-yM)^2 + (zE-zM)^2); %distance Earth-Moon


aE = -((Mu_Sun+Mu_Earth)/(pE*pE*pE)); %acceleration of the Earth due to the Sun
aMS = -((Mu_Sun+Mu_Moon)/(pM*pM*pM)); %acceleration of the Moon due to the Sun
aME = -((Mu_Earth+Mu_Moon)/(pEM*pEM*pEM)); %acceleration of the Moon due to the Earth

aSS = -((Mu_Sun)/(pS*pS*pS)); %acceleration of the satellite due to the Sun
aSE = -((Mu_Earth)/(pSE*pSE*pSE)); %acceleration of the satellite due to the Earth
aSM = -((Mu_Moon)/(pSM*pSM*pSM)); %acceleration of the satellite due to the Moon.



xaE = aE.* xE ; %accelaration in x-direction of the Earth
yaE = aE.* yE ; %acceleration in y-direction of the Earth
zaE = aE.* zE ; %acceleration in z-direction of the Earth
xaM = aMS* xM + aME*(xM-xE); %acceleration in x-direction of the Moon 
yaM = aMS* yM + aME*(yM-yE); %acceleration in y-direction of the Moon 
zaM = aMS* zM + aME*(zM-zE); %acceleration in z-direction of the Moon 
xaS = aSS.* xS + aSE*(xS-xE) + aSM*(xS-xM); %acceleration in  x-direction of the satellite
yaS = aSS.* yS + aSE*(yS-yE) + aSM*(yS-yM); %acceleration in y-direction of the satellite
zaS = aSS.* zS + aSE*(zS-zE) + aSM*(zS-zM); %acceleration in z-direction of the satellite

f = [xvE yvE zvE xaE yaE zaE xvM yvM zvM xaM yaM zaM xvS yvS zvS xaS yaS zaS]';