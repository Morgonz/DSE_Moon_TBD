function [ const ] = walkerdeltagenerator( i,n,p,h )
%walkerdeltagenerator: generate walker delta orbit parameters based on
%input requirements
const = [];
R_M = 1737.1e3;

% n = amount of orbit planes
% p = amount of sats per orbit
deltaRAAN = pi/n;
deltaphase = 2*pi/p;
e = 1;
omega = 0;
ID = 0;
for idx=1:n %per plane
    ID = idx*100;
   RAAN = deltaRAAN*(n-1);
    for jdx=1:p %per sat
        ID = ID + 1;
        phase = deltaphase*(p-1);
        
        append = [ID e i RAAN omega phase h]';
        const = [const append];
        keplerplot(append(2),append(3),append(4),append(5),append(6),append(7));
    end
end
[p,q,r] = sphere(40);
p = p.*R_M; q = q.*R_M; r = r.*R_M;

surf(p,q,r);hold on;
axis vis3d
axis equal
hold off;
        

end

