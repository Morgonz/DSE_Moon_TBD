nplanes=12;
A=zeros(7,12);
AAN0=0;
INC=60;
AP=80;
rp=(1736+700);
ra=9692;
DV=0;

%plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
for i=(1:nplanes)
    A(i,1)=AAN0+(360/nplanes)*i;
    A(i,2)=INC;
end

for i=(1:(nplanes-1))
    plane=plane_change(A(i,1),A(i,2),AP,ra,rp,A(i+1,2),A(i+1,2)).run();
    DV=DV+plane.dv;
end

DV