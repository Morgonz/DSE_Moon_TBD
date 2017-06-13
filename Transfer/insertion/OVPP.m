%One vehicle inserts in each plane, vehicles launched at the same time

nplanes=3;          %number of vehicles launched at once
rleo=6671;          %starting LEO radius
ra=10000;           %apoapse after initial manoeuvre
rp=1736+700;        %final moon orbit radius
vinf=0.830;
MuM=4.905*10^12;
Apl=zeros(nplanes,2);
for i=(1:nplanes)
    Apl(i,1)=(i-1)*360/12;
    Apl(i,2)=65;
end

%guess
ILON=25;
IINC=[62,75,80];

vp1=((vinf*10^3)^2+2*MuM/rp*10^(-3))^0.5;
vp2=(MuM*(2/rp*10^(-3)-2/(ra+rp)*10^(-3)))^0.5;
vc=(MuM/rp*10^(-3))^0.5;
dv1=vp1-vp2;
dv2=[0,0,0];
dv3=vp2-vc;

for i=(1:nplanes)
    %insertion_coords(W_int,phi_int,r_pea,v_inf)
    in=insertion_coords(ILON,IINC(i),rp,vinf).run();
    %plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
    plane=plane_change(in.W_an,in.inc,in.w_pea,ra,rp,Apl(i,1),Apl(i,2)).run();
    dv2(i)=plane.dv;
end

sum(dv2)
dv2
dv1+dv3
