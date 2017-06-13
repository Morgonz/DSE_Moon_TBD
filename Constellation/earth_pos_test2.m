    RAANe=1.4180;
    INCe=0.1177;
    theta_in=1.5835;
    an0=0.5988;
    e_m=0.0450;
    sma_m=3.8140e+08;
    
    an=tanomaly(an0,t,sma,e_m,.001);
    theta=an+theta_in;
    
    T3=rotationmatrices.Tz(RAANe);
    T2=rotationmatrices.Tx(INCe);
    T1=rotationmatrices.Tz(theta);
    
    E_coords=T3*T2*T1*[0;0;1]*sma_m*(1-e_m^2)/(1+e*cos(theta));