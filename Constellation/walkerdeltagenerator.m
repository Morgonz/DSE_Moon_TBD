function [  ] = walkerdeltagenerator(i,T,P,F,h)
%walkerdeltagenerator: generate walker delta orbit parameters based on

%input requirements:
% i = inclination [deg]
% n = amount of orbit planes
% p = amount of sats per orbit
% h = altitude from Moon surface

Rm = 1737.1e3;
e = 0;
omega = 0;
%i calculate from polar coverage
i = deg2rad(i);
const = [];
loc_sat = [];

n_sat = T/P;
%F = P-1-T/P; % F = 0:P-1
notice_F = [num2str(F) ' vs ' num2str(P-1-T/P)];
disp(notice_F)

% Orbit plane & sat phasing spacings
deltaRAAN = 2*pi/P
deltaphase = 2*pi/T*F

figure;
hold on;
for idx=1:P %per plane
    ID = idx*100;
   RAAN = deltaRAAN*(idx-1);
    for jdx=1:n_sat %per sat
        ID = ID + 1;
        phase = deltaphase*(jdx-1);
        
%         sat = [ID e i RAAN omega phase]';
%         const = [const sat];
        [X,Y,Z] = keplerplot(e,i,RAAN,omega,phase,h);
        scatter3(X(1),Y(1),Z(1),20,'filled');
        loc_sat = [loc_sat;X(1) Y(1) Z(1) ID];
    end
    plot3(X,Y,Z,'LineWidth',1);
end

% for k=1:(n_plane*p_sat)
%     RAAN = 2*pi*(k-1)/n_plane;
%     phase = 2*pi/(p_sat*n_plane)*(k-1)+2*pi/(p_sat)*(k-1)
% %    phase = 2*pi*(k-1)/n_plane;
%     [X,Y,Z] = keplerplot(e,i,RAAN,omega,phase,h);
%     scatter3(X(1),Y(1),Z(1),20,'filled');
%     plot3(X,Y,Z,'LineWidth',1.5);
%     loc_sat = [loc_sat;X(1) Y(1) Z(1)];
%     
% end


%plot moon sphere
[p,q,r] = sphere(40);
p = p.*Rm; q = q.*Rm; r = r.*Rm;

surf(p,q,r);hold on;
axis vis3d
axis equal
hold off;
% for k=1:nsat
%   % Right ascension of ascending node, in radians
%   thetao = ran0*deg2rad + 2*pi*(k-1)/nplanes;
% 
%   % True anomaly sequence, in radians
%   thetav = anom0*deg2rad + 2*pi*harmonic*(k-1)/nplanes + 2*pi*t/T;
% 
%   % Compute orbit from kepler elements
%   [XECI{k},VECI{k}] = kepl2eci(radius,0,thetai,thetao,0,thetav);
% end;

end

