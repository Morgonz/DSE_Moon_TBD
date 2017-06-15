function yrot = RotatingFrameSunEarth(state)
xE = state(:,1);
yE = state(:,2);
vxE = state(:,4);
vyE = state(:,5);
xS = state(:,7);
yS = state(:,8);
vxS = state(:,10);
vyS = state(:,11);

xEr =[];
yEr =[];
xSr =[];
ySr =[];
vxEr = [];
vyEr = [];
vxSr = [];
vySr = [];

rE=149.60E+9;
Mu_Earth = 3.98574405E+14; %m^3 s^-2
Mu_Sun = 1.327124E+20; %m^3 s^-2 from wikipedia

omega = 1/sqrt((rE^3)/(Mu_Earth+Mu_Sun));%rad/s

for i = linspace(1,length(xE),length(xE))
    xE1 = xE(i);
    yE1 = yE(i);
    if xE1>0
        tet = atan(yE1/xE1);
    elseif xE1<0
        tet = atan(yE1/xE1)+pi;
    else
        if yE1>0
            tet = pi/2;
        else
            tet = -pi/2;   
        end
    end
    xEr1 = xE(i)*cos(-tet)-yE(i)*sin(-tet);
    yEr1 = xE(i)*sin(-tet)+yE(i)*cos(-tet);
    xSr1 = xS(i)*cos(-tet)-yS(i)*sin(-tet);
    ySr1 = xS(i)*sin(-tet)+yS(i)*cos(-tet);
    vxEr1 = vxE(i)+(omega*yE(i));
    vyEr1 = vyE(i)-(omega*xE(i));
    vxSr1 = vxS(i)+(omega*yS(i));
    vySr1 = vyS(i)-(omega*xS(i));
    
    
    xEr = [xEr; xEr1];
    yEr = [yEr; yEr1];   
    xSr = [xSr; xSr1];
    ySr = [ySr; ySr1];
    vxEr = [vxEr; vxEr1];
    vyEr = [vyEr; vyEr1];
    vxSr = [vxSr; vxSr1];
    vySr = [vySr; vySr1];
end
yrot = [xEr yEr state(:,3) vxEr vyEr state(:,6) xSr ySr state(:,9) vxSr vySr state(:,12)];