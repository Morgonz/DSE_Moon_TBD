function yrot = RotatingFrame(state)
x = state(:,1);
y = state(:,2);
xM = state(:,7);
yM = state(:,8);

xr =[];
yr =[];
xMr=[];
yMr=[];
for i = linspace(1,length(x),length(x))
    xS = xM(i);
    yS = yM(i);
    if xS>0
        tet = atan(yS/xS);
    elseif xS<0
        tet = atan(yS/xS)+pi;
    else
        if yS>0
            tet = pi/2;
        else
            tet = -pi/2;   
        end
    end
    xr1 = x(i)*cos(-tet)-y(i)*sin(-tet);
    yr1 = x(i)*sin(-tet)+y(i)*cos(-tet);
    xMr1 = xS*cos(-tet)-yS*sin(-tet);
    yMr1 = xS*sin(-tet)+yS*cos(-tet);
    xr = [xr; xr1];
    yr = [yr; yr1];
    xMr = [xMr; xMr1];
    yMr = [yMr; yMr1];
end

yrot = [xr yr state(:,3) state(:,4) state(:,5) state(:,6) xMr yMr state(:,9) state(:,10) state(:,11) state(:,12)];
