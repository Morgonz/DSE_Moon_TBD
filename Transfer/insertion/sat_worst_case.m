nplanes=12;
A=zeros(12,2);
AAN0=90;
INC0=60;
INC=67;
AP=90;
rp=2353;
ra=9692;
MuM=4.905*10^12;

DVT=0;
DVPP=zeros(1,nplanes);

%plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
for i=(1:nplanes)
    A(i,1)=(360/nplanes)*i;
    A(i,2)=INC;
end

for i=(1:nplanes)
    a=plane_change(AAN0,INC0,AP,ra,rp,A(i,1),A(i,2)).run();
    DVT=DVT+a.dv;
    DVPP(i)=a.dv;
end

a1=(ra+rp)/2*10^3;
V_cir=(MuM/(rp*10^3))^(1/2);
DVC=(MuM*(2/(rp*10^3)-1/a1))^(1/2)-V_cir;

max(DVPP)
min(DVPP)
sum(DVPP)
