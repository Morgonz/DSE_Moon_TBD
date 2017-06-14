function dvt=TIM(aa)
    MuM=4.905*10^12;
    ra=(aa+1737+1629)*10^3;
    rp=(1737+1629)*10^3;
    r_f=(100+1737)*10^3;
    
    vc=(MuM/rp)^.5;
    v1p=(MuM*(2/rp-2/(ra+rp)))^.5;
    v1a=(MuM*(2/ra-2/(ra+rp)))^.5;
    v2a=(MuM*(2/ra-2/(ra+r_f)))^.5;
    
    dv1=v1p-vc;
    dv2=v1a-v2a;
    dvt=dv1+dv2;
end