function [f] = DeploymentSimulationAcc(t,state)

Mu_Moon = 4.902801e12; %m^3 s^-2
R_Moon = 1737400; %m mean Moon radius

xp = [];
yp = [];
zp = [];
vx = [];
vy = [];
vz = [];
for i = linspace(1,length(state)-5,(length(state)/6))
    xp = [xp; state(i)];
    yp = [yp; state(i+1)];
    zp = [zp; state(i+2)];
    vx = [vx; state(i+3)];
    vy = [vy; state(i+4)];
    vz = [vz; state(i+5)];
end

p = sqrt(xp.^2 + yp.^2 + zp.^2);
a = [];
for i = linspace(1,length(p),length(p))
    a1 = -(Mu_Moon/(p(i)^3));
    a  = [a; a1];
end

ax = a.*xp;
ay = a.*yp;
az = a.*zp;

f1 = [];
for i = linspace(1,length(a),length(a))
    f1 = [f1; vx(i); vy(i); vz(i); ax(i); ay(i); az(i)];
end
f = f1;

    