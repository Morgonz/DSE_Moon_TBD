%define (orbit)
nplanes=12;
MuM=4.905*10^12;
ra=()*10^3;
a2=(10000+1736+1629)/2*10^3;
a_eol=(1730*10^3+10000*10^3)/2;
INC=67*2*pi/360;
DRAAN=30*pi*2/360;
Di=360/nplanes*2*pi;
Vc=(MuM/(1736+700)*10^(-3))^0.5;
a1=(1736+1629)*10^3;

%define (propulsion)
Isp=3000;
F=10;
m1=10;

DVl=pi/2*(MuM/a)^0.5*abs(DRAAN)*sin(INC)*(nplanes-1);
DVi=Vc*(2-2*cos(pi/2*Di))^0.5;
DVa=(MuM/a2)^0.5-(MuM/a)^0.5;
DV_eol=-(MuM/a)^0.5+(MuM/a_eol)^0.5;

DTl=m1*Isp/F*(exp(-DVl)-1);