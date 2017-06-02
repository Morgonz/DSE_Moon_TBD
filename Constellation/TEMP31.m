%walkerdeltagenerator: generate walker delta orbit parameters based on
t = 1; % time form 1-100 for one orbit period
i = 67;
T = 72;
P = 6;
F = 1;
h = 1627e3;
%input requirements:
% i = inclination [deg]
% n = amount of orbit planes
% p = amount of sats per orbit
% h = altitude from Moon surface

Rm = 1738.1e3;
e = 0;
omega = 0;
%i calculate from polar coverage
i = deg2rad(i);
const = [];
loc_sat = [];
ID = 0;

n_sat = T/P;
%F = P-1-T/P; % F = 0:P-1
F_check = P-1-T/P;
notice_F = [num2str(F) ' vs ' num2str(F_check)];
disp(notice_F) 
% if F_check<0
%     F = F_check + P
% end
    

    

% Orbit plane & sat phasing spacings
deltaRAAN = 2*pi/P;
deltaphase = 2*pi/T*F;


hold on;
figure(1);
loc_sat = [];
for kdx=1:T
    RAAN = deltaRAAN*(kdx-1);
    phase = deltaphase*(kdx-1);

    [X,Y,Z] = keplerplot(e,i,RAAN,omega,phase,h);
    
    scatter3(X(t),Y(t),Z(t),25,'filled');
    plot3(X,Y,Z,'LineWidth',1.25);
    ID = ID + 1;
    loc_sat = [loc_sat;X(1) Y(1) Z(1) ID];

end


% for idx=1:P %per plane
%     ID = idx*100;
%    RAAN = deltaRAAN*(idx-1);
%     for jdx=1:n_sat %per sat
%         ID = ID + 1;
%         phase = deltaphase*((idx-1)*n_sat+(jdx-1));
%         
% %         sat = [ID e i RAAN omega phase]';
% %         const = [const sat];
%         [X,Y,Z] = keplerplot(e,i,RAAN,omega,phase,h);
%         scatter3(X(1),Y(1),Z(1),20,'filled');
%         loc_sat = [loc_sat;X(1) Y(1) Z(1) ID];
%     end
%     plot3(X,Y,Z,'LineWidth',1);
% end

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
set(gca,'FontSize',20)
set(gca,'XtickLabel',[],'YtickLabel',[]);
set(gca,'visible','off')

[p,q,r] = sphere(40);
p = p.*Rm; q = q.*Rm; r = r.*Rm;
xlabel('X [m]','FontSize',20);
ylabel('Y [m]','FontSize',20);
zlabel('Z [m]','FontSize',20);

surf(p,q,r);hold on;

box on
grid on
axis vis3d
axis equal
hold off;
view(-45,5);
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
