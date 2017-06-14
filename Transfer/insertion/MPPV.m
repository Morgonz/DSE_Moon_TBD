nplanes=6;          %number of planes to insert into
rleo=6671;          %starting LEO radius
ra=10000;           %apoapse after initial manoeuvre
rp=1736+1629;        %final moon orbit radius
reol=1720;
vinf=0.830;
MuM=4.905*10^12;

%constellation info
INCC=50.2;
AAND=360/6;
A=zeros(nplanes,2);
for i=(1:nplanes)
    A(i,1)=AAND*(i-1);
    A(i,2)=INCC;
end

%insertion_coords(W_int,phi_int,r_pea,v_inf)
%in=insertion_coords(0,-70,rp,vinf).run();
in=insertion_coords(A(1,1),A(1,2),rp,vinf).run();   %initial plane (at vinf) assumed to alligned with the first to insert into
w_p=in.w_pea;
%A(i,1)=in.W_an;
%A(i,2)=in.inc;

dv2=zeros(nplanes-1,1);

%plane changes
%plane_change(aanI,incI,aopI,apaI,peaI,aanF,incF)
for i=(1:(nplanes-1))
    plane=plane_change(A(i,1),A(i,2),w_p,ra,rp,A(i+1,1),A(i+1,2)).run();
    w_p=plane.aopF;
    dv2(i)=plane.dv;
end

vp1=((vinf*10^3)^2+2*MuM/rp*10^(-3))^0.5;
vp2=(MuM*(2/rp*10^(-3)-2/(ra+rp)*10^(-3)))^0.5;
vc=(MuM/rp*10^(-3))^0.5;
dv1=vp1-vp2;
dv3=vp2-vc;
va_eol=(MuM*(2/ra*10^(-3)-2/(ra+reol)*10^(-3)))^0.5;
va=(MuM*(2/ra*10^(-3)-2/(ra+rp)*10^(-3)))^0.5;
dv_eol=va-va_eol;

sum(dv2)
