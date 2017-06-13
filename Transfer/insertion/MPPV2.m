nplanes=6;          %number of planes to insert into
rleo=6671;          %starting LEO radius
ra=1737+1629;           %apoapse after initial manoeuvre
rp=1737+1629;        %final moon orbit radius
reol=1737+100;
vinf=0.830;
MuM=4.905*10^12;

Isp=330;
Msat=15.42;
Mbase=504;
m=4;
g=9.81;

%initial plane
INC0=52.5;
AAN0=-15;

%constellation info
INCC=52.5;
AAND=360/6;
A=zeros(nplanes,2);
for i=(1:nplanes)
    A(i,1)=AAND*(i-1);
    A(i,2)=INCC;
end

w_p=80;

dv2=zeros(nplanes-1,1);

%plane changes
%plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
%plane=plane_change(AAN0,INC0,w_p,ra,rp,A(1,1),A(1,2)).run();
%w_p=plane.aopF*360/2/pi;
for i=(1:(nplanes-1))
    plane=plane_change(A(i,1),A(i,2),w_p,ra,rp,A(i+1,1),A(i+1,2)).run();
    w_p=plane.aopF*360/2/pi;
    dv2(i)=plane.dv;
    plane.cplot()
end

%dvs
vp1=((vinf*10^3)^2+2*MuM/rp*10^(-3))^0.5;
vp2=(MuM*(2/rp*10^(-3)-2/(ra+rp)*10^(-3)))^0.5;
vc=(MuM/rp*10^(-3))^0.5;
dv1=vp1-vp2;
dv3=vp2-vc;
va_eol=(MuM*(2/ra*10^(-3)-2/(ra+reol)*10^(-3)))^0.5;
va=(MuM*(2/ra*10^(-3)-2/(ra+rp)*10^(-3)))^0.5;
dv_eol=va-va_eol;

%mass estimation
Mrat2=exp(dv2/Isp/g);
Mrat_eol=exp(dv_eol/Isp/g);
Mrat_ins=exp(dv1/Isp/g);
Mt=Mbase*Mrat_eol;
for i=(1:nplanes-1)
    Mt=(Mt+m*Msat)*Mrat2(i);
end
Mt=Mt*Mrat_ins;
C=[ra,Mt,dv3+dv_eol];
Mt