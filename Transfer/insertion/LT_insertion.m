%define (orbit)
MuM=4.905*10^12;
a1=(1736+1629)*10^3;
a2=(1736+700)*10^3;
incr1=65/180*pi;
incr2=27/180*pi;
a_eol=(1736+700+1730)/2*10^3;

Vc1=(MuM/a1)^0.5;
Vc2=(MuM/a2)^0.5;

%define (propulsion)
Isp=4190;
F=0.34;
m1=2400;

DV1=(Vc2^2+Vc2^2-2*Vc2*Vc2*cos(incr2*pi/2))^0.5;
DV2=(Vc1^2+Vc1^2-2*Vc1*Vc1*cos(incr2*pi/2))^0.5;
DV3=(Vc1^2+Vc2^2-2*Vc1*Vc2*cos(0))^0.5;
DV_eol=-(MuM/a2)^0.5+(MuM/a_eol)^0.5;

DT1=Isp/F*m1*(1-exp(-DV1/Isp))/60/60/24;
DT2=Isp/F*m1*(1-exp(-DV2/Isp))/60/60/24;
DT3=Isp/F*m1*(1-exp(-DV3/Isp))/60/60/24;

dv1=(MuM/a1)^0.5-(MuM*(2/a1-2/(a1+a2)))^0.5;
dv2=-(MuM/a2)^0.5+(MuM*(2/a2-2/(a1+a2)))^0.5;