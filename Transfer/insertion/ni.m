MuM=4.905*10^12;
MuE=3.986*10^14;
% r1=(1629+1736)*10^3;
% rp2=1735*10^3;
% 
% dv=(MuM/r1)^0.5-(MuM*(2/r1-2/(r1+rp2)))^0.5;

rm=60000*10^3;
re=384399*10^3-rm;

vesc=(2*MuM/rm+2*MuE/re)^0.5;