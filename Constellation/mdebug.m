%a=walker_maintenance(INC,n_planes,SMA,T);
tic
a=walker_maintenance(50.2*pi/180,6,(1629+1737)*10^3,60*60*24*7);
a.options=odeset('RelTol',1e-12);
a=a.prop(1);
toc

%%
plot(a.t{1},a.sma{1})
%%
plot(a.t{1},a.ecc{1})
