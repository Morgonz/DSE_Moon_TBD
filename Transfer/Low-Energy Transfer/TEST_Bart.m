clc;
clear all;

    xE = 0;
    yE = 1;
    xS=0;
    yS=2;
    
    omega=3;
    vxE=-omega;
    vyE=0;
    vxS=0;
    vyS=0;
    
    tet=atan2(yE,xE);
    
    xEr = xE*cos(-tet)-yE*sin(-tet);
    yEr = xE*sin(-tet)+yE*cos(-tet);
    xSr = xS*cos(-tet)-yS*sin(-tet);
    ySr = xS*sin(-tet)+yS*cos(-tet);
    
    vxEr = vxE*cos(-tet)-vyE*sin(-tet)+xE*omega*sin(-tet)+yE*omega*cos(-tet);
    vyEr = vxE*sin(-tet)+vyE*cos(-tet)-xE*omega*cos(-tet)+yE*omega*sin(-tet);
    
    vxSr = vxS*cos(-tet)-vyS*sin(-tet)+xS*omega*sin(-tet)+yS*omega*cos(-tet);
    vySr = vxS*sin(-tet)+vyS*cos(-tet)-xS*omega*cos(-tet)+yS*omega*sin(-tet);
