function [ loc_sat ] = walkerdeltagenerator(i,T,P,F,h,t)
%walkerdeltagenerator: generate walker delta orbit parameters based on
%input requirements:
% i = inclination [deg]
% T = amount of satellites
% p = amount of sats per orbit
% h = altitude from Moon surface [m]
% t = nondimensional time from 1-100

% i = 50.2;
% T = 36;
% P = 6;
% F = 1;
% h = 1627e3;
% t = 1;

Rm = 1738.1e3;
e = 0;
omega = 0;
%i calculate from polar coverage
i = deg2rad(i);
const = [];
loc_sat = [];
ID = 1;

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
    
    if ID==601||ID==602||ID==603||ID==604||ID==605
        ID = ID - 500 + 1;
        disp('akl')
    else
        ID = ID + 100;
    end

    loc_sat = [loc_sat;X(1) Y(1) Z(1) ID];

end

%plot moon sphere
set(gca,'FontSize',20)
set(gca,'XtickLabel',[],'YtickLabel',[]);
%set(gca,'visible','off')

[p,q,r] = sphere(40);
p = p.*Rm; q = q.*Rm; r = r.*Rm;
xlabel('X [m]','FontSize',12);
ylabel('Y [m]','FontSize',12);
zlabel('Z [m]','FontSize',12);

surf(p,q,r);hold on;

box on
grid on
axis vis3d
axis equal
hold off;
view(-45,5);


end