% This program is a test program to see the applications of the animation
% of a cannon ball

% clear all;
% close all;
% clc;
test =

%% Constants
% Moon Parameters
GMm = 4.902801e12;
Rm  = 1738100; % Equatorial radius (km) 1738.1

%% initial sat
alt_sat = 500e3;
[gx, gy, gz] = gravitysphericalharmonic([Rm+alt_sat 0 0],'LP165P',0);
GMm = abs(gx*(Rm+alt_sat)^2);
v_sat = sqrt(GMm/(Rm+alt_sat));
Ts = 2*pi*sqrt((Rm+alt_sat)^3/GMm); % sec
NT = 268;
T = NT*Ts;

%% Satellite parameters
x_sat = 0;
y_sat = 0;
z_sat = Rm + alt_sat;

vx_sat = 0;
vy_sat = v_sat;
vz_sat = 0;

%% Integration
% options of the integrator
options = odeset('RelTol',1e-6);

% The time integration of the sat
[t0,y0] = ode45(@sat2BPdeg0,[0 T],[x_sat y_sat z_sat vx_sat vy_sat vz_sat],options);
[tH,yH] = ode45(@sat2BPdegH,[0 T],[x_sat y_sat z_sat vx_sat vy_sat vz_sat],options);

[Xm,Ym,Zm] = sphere(40);  % Equatorial radius (km) 1738.1

%% Plotting of orbits
disc_yH = round(length(yH)/NT);

figure
surf(Xm*Rm, Ym*Rm, Zm*Rm)
hold on
plot3(y0(:,1),y0(:,2),y0(:,3),'g')
plot3(yH(:,1),yH(:,2),yH(:,3),'r')
plot3(yH(end-disc_yH:end,1),yH(end-disc_yH:end,2),yH(end-disc_yH:end,3),'p')
%scatter3(3000e3,0,0,'p')
%scatter3(0,3000e3,0, 'y')
title('orbit evolution 3d')

axis vis3d
axis equal

r0 = sqrt(y0(:,1).^2 + y0(:,2).^2 + y0(:,3).^2);
rH = sqrt(yH(:,1).^2 + yH(:,2).^2 + yH(:,3).^2);

%% Radius evolution of orbit
figure
plot(t0,r0,'g')
hold on
plot(tH,rH,'r')
plot([0,tH(end)],[Rm,Rm])
title('radius evolution')
hold on
myfit = polyfit(tH, rH-rH(1), 1);
plot(tH, rH(1)+tH*myfit(1)+myfit(2));
maxdiff_rH = max(rH - mean(r0))

%% Evolution of semi-major axis
% V^2 = µ(2/r - 1/a)
V0 = sqrt(y0(:,4).^2 + y0(:,5).^2 + y0(:,6).^2);
VH = sqrt(yH(:,4).^2 + yH(:,5).^2 + yH(:,6).^2);

a0 = -1./((V0.^2)/GMm - 2./r0);
aH = -1./((VH.^2)/GMm - 2./rH);

figure
plot(t0(1:end-1),a0(1:end-1)-a0(1),'g')
hold on
plot(tH(1:end-1),aH(1:end-1)-a0(1),'r')
title('semi-major axis evolution')
hold on
myfit = polyfit(tH(1:end-1), aH(1:end-1),1);
plot(tH, myfit(2)+tH*myfit(1)-a0(1));

maxdiff_aH = max(aH - mean(a0))

%% Evolution of eccentricity
% e = (ra-rp)/(ra+rp)
rp = [];
idp = [];
ra = [];
vp = [];
ida = [];
va = [];
te = [];
disc_rH = length(rH)/NT;
for idx=1:NT
    rp = [rp min(rH(1+round((idx-1)*disc_rH):round(idx*disc_rH)))];
    idp = [idp max(find(rp(end)==rH))];
    ra = [ra max(rH(1+round((idx-1)*disc_rH):round(idx*disc_rH)))];
    ida = [ida max(find(rp(end)==rH))];
    adx = find(rH==rp(end));
    pdx = find(rH==ra(end));
    vp = [vp VH(pdx(end))];
    va = [va VH(adx(end))];
    te = [te idx];       
end

e = (ra-rp)./(ra+rp);

figure
plot(te,ra,'g')
hold on
plot(te,rp,'r')
title('ra (g) and rp (r) evolution')

figure
plot(te,e)
title('eccentricity evolution')

%% Inclination
u_a = yH(ida,1:3);
u_p = yH(idp,1:3);

u_az = u_a;
u_pz = u_p;

u_az(:,3) = 0;
u_pz(:,3) = 0;

incl = [];
for idx=1:length(u_a)
    ang = atan2(norm(cross(u_a(idx,:),u_az(idx,:))) , dot(u_a(idx,:),u_az(idx,:)));
    incl = [incl rad2deg(ang)];
end

figure
plot(te,incl)
title('inclination evolution')

%% Delta v calculation
% https://www.google.nl/imgres?imgurl=http%3A%2F%2Fwww.braeunig.us%2Fspace%2Fpics%2Feq4-58.gif&imgrefurl=http%3A%2F%2Fwww.braeunig.us%2Fspace%2Forbmech.htm&docid=bu4nm38VdM1DrM&tbnid=robluQbRmBPU_M%3A&vet=10ahUKEwjtq_Se5vnTAhWIbVAKHdn6AGIQMwgnKAIwAg..i&w=555&h=313&bih=686&biw=1536&q=energy%20equation%20orbit%20equation&ved=0ahUKEwjtq_Se5vnTAhWIbVAKHdn6AGIQMwgnKAIwAg&iact=mrc&uact=8

r_a = rp(end);
r_p = Rm+alt_sat;
a_tx = (r_a + r_p)/2;

v_ia = va(end);
%v_ia = sqrt(GMm/r_a);
v_fb = v_sat;

v_txa = sqrt(GMm*(2/r_a - 1/a_tx));
v_txb = sqrt(GMm*(2/r_p - 1/a_tx));

d_va = v_ia - v_txa;
d_vb = v_fb - v_txb;
d_vt = abs(d_va) + abs(d_vb)

maint_s = d_vt/T;
main_y  = maint_s*365*24*3600




