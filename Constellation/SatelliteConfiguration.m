%% Initial parameters
%Constellation inputs
%RP = 2; %Redundancy parameter [-]
h = 700e3; %Network satellite height [m]
i = deg2rad(10); %Elevation angle

S_FOV = deg2rad(120/2); %Satellite field of view
R_FOV = deg2rad(120/2); %Rover field of view

%Moon parameters5
%M_moon = 7.342e22; % [kg]
%G = 6.67408e-11;
Rm = 1737.1e3; % [m]

Rs = Rm + h;


%Ground swath width
a = asin(Rm/Rs); %half angle of satellite vision on the Moon
%b = acos(R_M*cos(i)/(R_S))-i; %half angle 

if a>S_FOV
    %Calculate moon coverage angle with ANTENNA FOV
    %b=deg2rad(90)-acos((R_S*sin(deg2rad(S_FOV)))/R_M);
    b = deg2rad(90) - S_FOV - acos((Rs*sin(S_FOV))/Rm);
    disp('Antenna limited FOV')

else
    %Calculate moon coverage angle with MOON FOV + elevation angle (moon surface roughness)
    b = acos(Rm*cos(i)/(Rs))-i;
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