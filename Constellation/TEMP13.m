e = 0;
RAAN = 0;
omega = 0;
dphase = 2*pi/3;
h = 6500e3;
i = [52 82];
loc_sat = [];
t = 1;
for idx=1:3
    phase = (idx-1)*dphase;
    [X,Y,Z] = keplerplot(e,i(1),RAAN,omega,phase,h);
    loc_sat = [loc_sat;X(t) Y(t) Z(t)];
    scatter3(X(t),Y(t),Z(t),25,'filled');hold on;
    plot3(X,Y,Z,'LineWidth',1.25);hold on;
    phase = (idx-1)*dphase+0.25*dphase;
    [X,Y,Z] = keplerplot(e,i(2),RAAN,omega,phase,h);
    loc_sat = [loc_sat;X(t) Y(t) Z(t)];
    scatter3(X(t),Y(t),Z(t),25,'filled');hold on;
    plot3(X,Y,Z,'LineWidth',1.25);hold on;
end




% set(gca,'FontSize',20)
% [p,q,r] = sphere(40);
% p = p.*Rm; q = q.*Rm; r = r.*Rm;
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% zlabel('Z [m]','FontSize',20);
% for i=1:length(loc_sat)
%     
%     scatter3(X(t),Y(t),Z(t),25,'filled');
%     plot3(X,Y,Z,'LineWidth',1.25);
% end
% 
% surf(p,q,r);hold on;
% box on
% grid on
% axis vis3d
% axis equal
% hold off;
% view(-45,5);