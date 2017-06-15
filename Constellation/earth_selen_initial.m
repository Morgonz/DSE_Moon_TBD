%% Ephemeris model and time
ep_model='432t';
ep_time0=2459001.500000000;
t=0;
ep_time=ep_time0+t/(60*60*24);

%% Initial Moon position and vel (ICRF, earth mean equator)
% ecliptic:
% rmx=-3.638790674623008e8; rmy=-3.401444133974022e7; rmz=3.338839823502487e7;
% vmx=1.243026318251819e2; vmy=-1.060865683273835e3; vmz=-1.317638329677234e1;
% mean eq:
% rmx=-3.638790674623008e8; rmy=-4.448878187156640e7; rmz=1.710308872512287e07;
% vmx=1.243026318251819e2; vmy=-9.680839703952510e2; vmz=-4.340772296361691e2;

[pos,vel]=planetEphemeris(ep_time,'Moon','Earth',ep_model,'km');
rmx=pos(1)*10^3; rmy=pos(2)*10^3; rmz=pos(3)*10^3;
vmx=vel(1)*10^3; vmy=vel(2)*10^3; vmz=vel(3)*10^3;

kepler=getkepler_e(rmx,rmy,rmz,vmx,vmy,vmz);
inc0=acos(kepler.nor(3));

%% Rotation matrix: ICRF(mean earth eq.) => Selenocentric
angles=moonLibration(ep_time,ep_model);
phi_lib=angles(1);
theta_lib=angles(2);
psi_lib=angles(3);
T_sel=rotationmatrices.Tz(psi_lib)*rotationmatrices.Tx(theta_lib)*rotationmatrices.Tz(phi_lib);

%% Earth position wrt Moon (selenocentric)
Re=[-rmx;-rmy;-rmz];
Re=T_sel*Re;
mRe=(Re(1)^2+Re(2)^2+Re(3)^2)^.5;

%% E-M plane normal vector
n_EM_ICRF=kepler.nor;
n_EM_sel=T_sel*n_EM_ICRF';

%% Ascending node vector
anvec=cross([0;0;1],n_EM_sel);

%% Inclination and long of scending node
inc=acos(n_EM_sel(3));
raan=acos(anvec(1)/(anvec(1)^2+anvec(2)^2+anvec(3)^2)^.5);
if anvec(2)<0
    raan=2*pi-raan;
end

%% Initial Earth position theta
theta_sel_in=acos((Re(1)*anvec(1)+Re(2)*anvec(2)+Re(3)*anvec(3))/mRe);
if Re(2)<0
    theta_sel_in=2*pi-theta_sel_in;
end

%% Initial Moon orbit true anomaly
an_in=kepler.ta;
