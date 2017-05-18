function [  ] = walkerdeltagenerator( n,p,h )
%walkerdeltagenerator: generate walker delta orbit parameters based on

%input requirements:
% i = inclination [deg]
% n = amount of orbit planes
% p = amount of sats per orbit
% h = altitude from Moon surface

Rm = 1737.1e3;
e = 0;
omega = 0;

const = [];
loc_sat = [];


% Orbit plane & sat phasing spacings
deltaRAAN = 2*pi/n;
deltaphase = 2*pi/p;

i = p/deltaRAAN
figure;
hold on;
for idx=1:n %per plane
    ID = idx*100;
   RAAN = deltaRAAN*(idx-1);
    for jdx=1:p %per sat
        ID = ID + 1;
        phase = deltaphase*(jdx-1);
        
%         sat = [ID e i RAAN omega phase]';
%         const = [const sat];
        [X,Y,Z] = keplerplot(e,i,RAAN,omega,phase,h);
        scatter3(X(1),Y(1),Z(1),20,'filled');
        loc_sat = [loc_sat;X(1) Y(1) Z(1) ID];
    end
    plot3(X,Y,Z,'LineWidth',1.5);
end
%plot moon sphere
[p,q,r] = sphere(40);
p = p.*Rm; q = q.*Rm; r = r.*Rm;

surf(p,q,r);hold on;
axis vis3d
axis equal
hold off;


end

