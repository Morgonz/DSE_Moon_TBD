function [f] = EarthAcc(t,state)
x=state(1);
y=state(2);
z=state(3);
vx=state(4);
vy=state(5);
vz=state(6);

Mu_Earth=3.9860E+14;
Mu_Sun=1.32712E+20;

r=sqrt(x^2 + y^2 + z^2);
a=-((Mu_Earth+Mu_Sun)/r^3);

ax=a*x;
ay=a*y;
az=a*z;

f=[vx vy vz ax ay az]';


