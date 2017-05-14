%% Initial parameters
%Constellation inputs
%R = 3; %Redundancy parameter [-]
h = 400e3; %Network satellite height [m]
i = deg2rad(5); %Elevation angle

S_FOV = deg2rad(120/2); %Satellite field of view
R_FOV = deg2rad(120/2); %Rover field of view

%Moon parameters
%M_moon = 7.342e22; % [kg]
%G = 6.67408e-11;
R_M = 1737.1e3; % [m]

R_S = R_M + h;



%Ground swath width
a = asin(R_M*cos(i)/(R_S));
b = acos(R_M*cos(i)/(R_S))-i;

if a>S_FOV
    %Calculate moon coverage angle with antenna FOV
    b=deg2rad(30)-acos((R_S*sin(deg2rad(60)))/R_M);
else
    %Calculate moon coverage angle with Moon FOV + elevation angle (moon surface roughness)
    b = acos(R_M*cos(i)/(R_S))-i;
end

%Ascending node orbit
dn = b; %Ascending node spacing
n = ceil(pi/dn); %Number of orbit planes
%phasing
dp = b; %Phase spacing
p = ceil(2*pi/dp); %Number of sats per plane

%satellite list
Nsats = p*n;
out = [num2str(n),' number of orbit planes, ',num2str(p),' sats per orbit, ',num2str(Nsats),' total satellites'];
disp(out)


%% Redundancy check
