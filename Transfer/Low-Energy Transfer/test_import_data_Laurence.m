y_inertial_bart = importdata('State_Earth_to_L2_PILSBAAS_LAURENCE.mat')
data_Bart = importdata('TIME_Earth_to_L2_backwards_Bart.mat')
yrot_bart = RotatingFrameSunEarth(y_inertial_bart);

hold on
plot3(yrot_bart(:,7),yrot_bart(:,8),yrot_bart(:,9),'r-','DisplayName','halo orbit');


Earth = plot3(rE,0,0,'bo','DisplayName','Earth location')
L2 = plot3(rE+rL,0,0,'k*','DisplayName','Sun-Earth L2')
%surf(X*R_Sun,Y*R_Sun, Z*R_Sun,'DisplayName','Sun position')
hold off
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
%legend([L2 Earth],{'Earth-Sun L2','Earth'})
legend('show')
axis equal
axis vis3d