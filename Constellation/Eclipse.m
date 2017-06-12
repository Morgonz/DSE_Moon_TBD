function [ T,T_eclipse,cycles ] = Eclipse( h,lifetime )
%Eclipse Calculates orbit time, eclipse time and cycles around the Moon
%   Detailed explanation goes here
%Moon parameters
M_moon = 7.342e22; % [kg]
G = 6.67408e-11;
Rm = 1738.1e3; % [m]
% d_SE = 149597871e3; % [m]

GMm = G*M_moon;
Rs = Rm + h; % [m]

%eclipse time
T = 2*pi*sqrt(Rs^3/GMm);
eclipse_angle = 2*(asin(Rm/Rs));
T_eclipse = eclipse_angle*sqrt(Rs^3/GMm);
%T_eclipse = 0; %for halo sats

cycles = (lifetime*365*24*3600)/T;

end

