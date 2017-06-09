function [dv1, dv2, dvt]=Hoh2(r0, inc0, rf)
% CR3B EMTransfer arriving at node

% define prking orbit
%   r0 [km] 
%   inc [deg] (ecliptic)
% define final orbit
%   rf [km]
r0=r0*10^3;
rf=rf*10^3;
inc0=inc0*2*pi/360;

%orbit parameters
MuE=3.986*10^14;
MuM=4.905*10^12;
incm=5.145/360*2*pi;
inc=inc0-incm;
rm=384399*10^3;
vm=(MuE/rm)^(1/2);
aT=(r0+rm)/2;
%velocities
v0=(MuE/r0)^(1/2);
vT1=(MuE*(2/r0-1/aT))^(1/2);
vT2=(MuE*(2/rm-1/aT))^(1/2);
vinf=((vm-vT2*cos(inc))^2+(vT2*sin(inc))^2)^(1/2);
vf1=(vinf^2+2*MuM/rf)^(1/2);
vf2=(MuM/rf)^(1/2);
%DVs
dv1=vT1-v0;
dv2=vf1-vf2;
dvt=dv1+dv2;
end

%[dv1, dv2, dvt]=Hoh2(6671,5.145,1837.1);