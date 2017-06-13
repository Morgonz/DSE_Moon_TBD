n=50;
Earth_coords=zeros(n,3);
r0_test=zeros(n,3);
v0_test=zeros(n,3);
t1=60*60*24*27;
tl=zeros(n,1);
sma=zeros(n,1);
ecc=zeros(n,1);
inc=zeros(n,1);
raan=zeros(n,1);
rx=zeros(n,1); ry=zeros(n,1); rz=zeros(n,1);
vx=zeros(n,1); vy=zeros(n,1); vz=zeros(n,1);
for i=(1:n)
    t=i/n*t1;
    tl(i)=t;
    
    % Ephemeris model and time
    ep_model='432t';
    ep_time0=juliandate(1990,12,1);
    ep_time=ep_time0+t/(60*60*24);

    % Rotation matrix: ICRF => Selenocentric
    angles=moonLibration(ep_time,ep_model);
    phi=angles(1);
    theta=angles(2);
    psi=angles(3);
    T_sel=rotationmatrices.Tz(phi)*rotationmatrices.Tx(theta)*rotationmatrices.Tz(psi);

    % Rotation matrix: Selenocentric => Inertial
    ome = 2*pi/(29.5*3600*24); %rad/s - https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
    OME  = -ome*t;
    T_I = rotationmatrices.Tz(OME);

    
    [pos,vel]=planetEphemeris(ep_time,'Moon','Earth','405','km');
    
    pos=T_I*T_sel*pos'; vel=T_I*T_sel*vel';
    
    rx(i)=pos(1)*10^3; ry(i)=pos(2)*10^3; rz(i)=pos(3)*10^3;
    vx(i)=vel(1)*10^3; vy(i)=vel(2)*10^3; vz(i)=vel(3)*10^3;
    r0_test(i,:)=[rx(i);ry(i);rz(i)];
    v0_test(i,:)=[vx(i);vy(i);vz(i)];
    R0_test(i)=(rx(i)^2+ry(i)^2+rz(i)^2)^.5;
    V0_test(i)=(vx(i)^2+vy(i)^2+vz(i)^2)^.5;
end

kepler=getkepler_e(rx,ry,rz,vx,vy,vz);

%% plot
%scatter3(Earth_coords(:,1),Earth_coords(:,2),Earth_coords(:,3))
scatter3(r0_test(:,1),r0_test(:,2),r0_test(:,3),'x')
axis vis3d;
axis equal;

%%
scatter(tl,kepler.SMA,'x')

%%
scatter(tl,R0_test,'x')

%%
scatter(tl,V0_test,'x')

%%
scatter(tl,kepler.ECC,'x')

%%
scatter(tl,rad2deg(kepler.INC),'x')

%%
scatter(tl,kepler.RAAN,'x')

%%
scatter(tl,kepler.AOP,'x')

%%
scatter(tl,kepler.nor(:,1),'x')
