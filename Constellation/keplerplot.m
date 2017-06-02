function [X,Y,Z] = keplerplot(e,inc,RAAN,omega,phase,h)
%keplerplot: Translate kepler input into 3d plot [X,Y,Z]

Rm = 1738.1e3; %[m]

Rs = Rm+h; %[m] minimum altitude (periapsis)
a = Rs/(1-e);

%from kep2chart_2002

theta=(linspace(phase,2*pi+phase,100))';
r = a*(1-e^2)./(1+e.*cos(theta));

%Coordinate generation
X = r.*(cos(RAAN).*cos(omega+theta)-sin(RAAN).*sin(omega+theta)*cos(inc));
Y = r.*(sin(RAAN).*cos(omega+theta)+cos(RAAN).*sin(omega+theta)*cos(inc));
Z = r.*(sin(inc).*sin(omega+theta));

% %% Uncomment for single orbit plot
% [p,q,r] = sphere(40);
% p = p.*Rm; q = q.*Rm; r = r.*Rm;
% 
% surf(p,q,r);hold on;
% plot3(X,Y,Z,'LineWidth',1.5);hold on;
% scatter3(X(1),Y(1),Z(1),16);hold on;
% axis vis3d
% axis equal

end

