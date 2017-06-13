%plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
%a=plane_change(0,60,-90,10000,700+1736,30,60).run().cplot();
a=plane_change(0,67,-90,10000,700+1736,30,67).run();
%a.aopF=-a.aopF2
a.cplot();

a.aopI*360/2/pi
a.aopF*360/2/pi