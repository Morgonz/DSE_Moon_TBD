%% Ephemeris model and time
ep_model='405';
ep_time0=juliandate(2030,1,1,0,0,0);
ep_time=ep_time0+t/(60*60*24);

%% Rotation matrix: ICRF => Selenocentric
angles=moonLibration(ep_time,ep_model);
phi=angles(1);
theta=angles(2);
psi=angles(3);
T_sel=rotationmatrices.Tz(phi)*rotationmatrices.Tx(theta)*rotationmatrices.Tz(psi);

%% E-M plane normal vector
n_EM_ICRF=[-0.0805;-0.4360;0.8963];
n_EM_sel=T_sel*n_EM_ICRF;

%% Ascending node vector
anvec=cross([0;0;1],n_EM_sel);

%% Inclination and long of scending node
inc=acos(n_EM_sel(3)/(n_EM_sel(1)^2+n_EM_sel(2)^2+n_EM_sel(3)^2)^.5);
raan=acos(anvec(1)/(anvec(1)^2+anvec(2)^2+anvec(3)^2)^.5);
