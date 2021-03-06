Mu_Earth=3.9860E+14;
Mu_Sun=1.32712E+20;
R_Sun=695700000;
R_Earth=6371000;

rE=1.496E+11;
vE=sqrt((Mu_Earth+Mu_Sun)/rE);

options=odeset('RelTol',2.22045e-14,'Events',@Cross);

[t,y,te]=ode113(@EarthAcc,[0 10^10],[rE 0 0 0 vE 0],options);
plot(y(:,1),y(:,2))
axis equal

function [value,isterminal,direction] = Cross(t,y)
value=y(2);
isterminal = 0;
direction = 1;
end