e = 0;
RAAN = 0;
omega = 0;
dphase = 2*pi/3;
h = 6500e3;
i = [52 82];
loc_sat = [];
t = 1;
for idx=1:3
    phase = (idx-1)*dphase;
    [X,Y,Z] = keplerplot(e,i(1),RAAN,omega,phase,h);
    loc_sat = [loc_sat;X(t) Y(t) Z(t)];
    phase = (idx-1)*dphase+0.25*dphase;
    [X,Y,Z] = keplerplot(e,i(2),RAAN,omega,phase,h);
    loc_sat = [loc_sat;X(t) Y(t) Z(t)];
end