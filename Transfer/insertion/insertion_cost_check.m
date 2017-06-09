rp=(1736+700)*10^3;
ra=10000*10^3;
vinf=0.8*10^3;
MuM=4.905*10^12;
inc=30*2*pi/360;

vp1=(vinf^2+2*MuM/rp)^0.5;
vc=(MuM/rp)^0.5;
DVT0=((vc*sin(inc))^2+(vp1-vc*cos(inc))^2)^0.5;

vp2=(MuM*(2/rp-2/(ra+rp)))^0.5;
va2=(MuM*(2/ra-2/(ra+rp)))^0.5;
DV1=vp1-vp2;
DV2=2*sin(inc/2)*va2;
DV3=vp2-vc;
DVT1=DV1+DV2+DV3;

DVT1/DVT0
