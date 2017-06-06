function yrot = RotatingFrameSunEarthMoon(state)
xE = state(:,1);
yE = state(:,2);
xM = state(:,7);
yM = state(:,8);
xS = state(:,13);
yS = state(:,14);

xEr =[];
yEr =[];
xMr =[];
yMr =[];
xSr =[];
ySr =[];

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
    xMr1 = xM(i)*cos(-tet)-yM(i)*sin(-tet);
    yMr1 = xM(i)*sin(-tet)+yM(i)*cos(-tet);
    xSr1 = xS(i)*cos(-tet)-yS(i)*sin(-tet);
    ySr1 = xS(i)*sin(-tet)+yS(i)*cos(-tet);
    xEr = [xEr; xEr1];
    yEr = [yEr; yEr1];
    xMr = [xMr; xMr1];
    yMr = [yMr; yMr1];    
    xSr = [xSr; xSr1];
    ySr = [ySr; ySr1];
end

yrot = [xEr yEr state(:,3) state(:,4) state(:,5) state(:,6) xMr yMr state(:,9) state(:,10) state(:,11) state(:,12) xSr ySr state(:,15) state(:,16) state(:,17) state(:,18)];
