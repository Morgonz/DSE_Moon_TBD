%define
n1=[0,1,1];
n2=[0,1,0];
SMA=1900;
SMA=SMA*10^3;
e=0.2;
%
n3=cross(n1,n2);
TA1=asin((n3(1)^2+n3(3)^2)^(1/2)/abs(n3(2)));
TA2=asin((n3(1)^2+n3(3)^2)^(1/2)/abs(n3(2)))+pi/2;
if TA1>=pi
    TA1=TA1-pi;
end
if TA2>=pi
    TA2=TA2-pi;
end
rm=SMA*(1-e^2)/(1-e*cos(TA2));


