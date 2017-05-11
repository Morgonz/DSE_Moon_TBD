% Sara Powell, ASEN 5050 Project. 

% This script calculates the equations of motions  for the third body
% according to Parkers book: Low energy Lunar Trajectory Design, 2014. 

function [result] = CRTBPmodelNOSTM_deltav(times,r,mu)
% Change variable names for understanding of the code. 
% Positions
x = r(1);
y = r(2);
z = r(3);
% Velocities
dotx = r(4);
doty = r(5);
dotz = r(6);

% Distance from the third body to the larger and smaller body respectively. Parker, 2014 page 36.
r1 = sqrt((x+mu)^2 + y^2 + z^2); 
r2 = sqrt((x-1+mu)^2 + y^2 + z^2); 

% Accelerations 
dotdotx = 2*doty + x -(1-mu)*((x+mu)/(r1^3)) - mu*(x-1+mu)/(r2^3);
dotdoty = -2*dotx + y - (1-mu)*(y/(r1^3)) - mu*(y)/(r2^3);
dotdotz = -(1-mu)*(z)/(r1^3) - mu*(z)/(r2^3); 

% create results vector. 
result = [r(4:6); [dotdotx ;dotdoty; dotdotz]];

end