% Sara Powell, ASEN5050 Project

% This code aims to develop a transfer from Earth to a halo Orbit
% around the L1 point. 


clear all;close all;clc
format long
addpath('D:\DSE\Functions')
% Standard parameters.
rad = (pi/180);
deg = (180/pi);
R_E = 6378136.3; %Earth Radius in meters
m_Earth = 5.9742e24; %kg
m_Sun = 1.9891e30;% kg
m_Moon = 7.34767309e22;% kg
% mu = (m_Earth+m_Moon)/(m_Sun + (m_Earth+m_Moon)); % Make sure this is the Moon-Earth barycenter!
mu = m_Moon/(m_Earth + m_Moon);
EM_dist = 384400;
parkingorbit_Earth = 185/EM_dist; %m LEO parking orbit. 

%% Plot initial orbit
% Initialize the state transition matrix
STM = reshape(eye(6),1,[]);
X_state = [0.8122;0;0;0;0.248312;0;STM'];
% Period of the orbit
orbitPeriod = 2.927456032136278;

% Propagate the orbit over time. 
ode_options = odeset('RelTol',1e-13,'AbsTol',1e-22);
[t_prop,X_prop] = ode45(@CRTBPmodel, [0 orbitPeriod], X_state, ode_options,mu);

figure(1)
plot3((X_prop(:,1)),X_prop(:,2),X_prop(:,3));hold on
plot3(0-mu,0,0,'ok', 'markerfacecolor', 'b', 'markersize', 22);hold on
plot3(1-mu,0,0,'ok', 'markerfacecolor', 'y', 'markersize',10);hold on
plot3(0.8369151324,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 8);hold on
plot3(1.1556821603,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 8);hold on
title('Initial Halo Orbit');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
grid on
grid on

% Just plot the orbit. 
figure(3)
plot3((X_prop(:,1)),X_prop(:,2),X_prop(:,3))
title('Halo Orbit','interpreter','latex','fontsize',30);
grid on;


%% Now use the single shooting differential correction
% Page 68 in Parker's book

% Single shooting method for constructing simple periodic symmetrix orbits
% in the CRTBP. (Used to construct halo orbits. )

% Orbit is symmetric about the y=0 plane and pierce the y=0 plane twice per
% orbit. X0 is the state of a simple periodic symmetric orbit at the y=0
% crossing with a positiive velocity in y(doty). and Thalf to be the state
% of the orbit at the one half orbit later when the y-velocity is negative
% at the crossing. 

% Change the options to find y=0 and stop integration at that point. 
ode_options1 = odeset('Events',@Findyzero,'RelTol',1e-13,'AbsTol',1e-13);
X0 = X_state;
% propagate the orbit until the half orbit. (Since y=0 at start, the next
% time y=0 means we have propagated half an orbit.)
deltas = [100;100];
cnt = 0;
while abs(norm(deltas))>1e-13
    cnt = cnt+1;
    [t_half,X_half] = ode45(@CRTBPmodel, [0 Inf], X0, ode_options1,mu);

    % Get the propagated state at T(1/2)
    dX = X_half(end,1:6);
    % Get the STM for this time.
    dPhi = reshape(X_half(end,7:end),6,[]);

    %% Compute the time derivative of the state. dX/dt
    dotx = dX(4);
    doty = dX(5);
    dotz = dX(6);

    % Distance from the third body to the larger and smaller body respectively. Parker, 2014 page 36.
    r1 = sqrt((dX(1)+mu)^2 + dX(2)^2 + dX(3)^2); 
    r2 = sqrt((dX(1)-1+mu)^2 + dX(2)^2 + dX(3)^2); 

    % Accelerations 
    dotdotx = 2*doty + dX(1) -(1-mu)*((dX(1)+mu)/(r1^3)) - mu*(dX(1)-1+mu)/(r2^3);
    dotdoty = -2*dotx + dX(2) - (1-mu)*(dX(2)/(r1^3)) - mu*(dX(2))/(r2^3);
    dotdotz = -(1-mu)*(dX(3))/(r1^3) - mu*(dX(3))/(r2^3); 
    % time derivative of the propagated state.
    dXdt = [dotx doty dotz dotdotx dotdoty dotdotz];
    %% Find the desired change in the final state's components should be. 
    % Only desirable to change dotx and dotz, but not important if the other 
    % components of the final state change. However, y will always be 0 since 
    % we propagate the orbit to T(1/2)

    CrazySTM = [dPhi(4,3)-dPhi(2,3)*(dotdotx/doty),dPhi(4,5)-dPhi(2,5)*(dotdotx/doty);...
        dPhi(6,3)-dPhi(2,3)*(dotdotz/doty), dPhi(6,5)-dPhi(2,5)*(dotdotz/doty)];

    deltas = CrazySTM^(-1)*[-X_half(end,4);-X_half(end,6)];
    norm(deltas)

    dXThalf = [0,0,deltas(1),0,deltas(2),0]';

    X0(1:6) = X0(1:6)+dXThalf;
end

%% Propagate and plot the new stable halo orbit. 
% Propagate the orbit over time. 
ode_options = odeset('RelTol',1e-13,'AbsTol',1e-22);
[t_halo,X_halo] = ode45(@CRTBPmodel, [0 t_half(end)*2], X0, ode_options,mu); 

% Plot the orbit togethet with the points. 
figure(4)
plot3((X_halo(:,1)),X_halo(:,2),X_halo(:,3),'linewidth',1.5);hold on
plot3(0-mu,0,0,'ok', 'markerfacecolor', 'b', 'markersize', 22);hold on
text(0-mu,0,0,'     Earth')
plot3(1-mu,0,0,'ok', 'markerfacecolor', 'y', 'markersize',10);hold on
plot3(0.8369151324,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 8);hold on
plot3(1.1556821603,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 8);hold on
title('Periodic Halo Orbit');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
grid on
grid on

figure(5)
plot3((X_halo(:,1))*EM_dist,X_halo(:,2)*EM_dist,X_halo(:,3)*EM_dist)
title('Halo Orbit','interpreter','latex','fontsize',30);
grid on;

%% Methodology page 125, Parker.
% 1) Construct the desired Halo Orbit. 
% 2) Construct the manifold segment
%     - Choose a tau valuem that is a point along the halo orbit. Choose a direction
%        and choose Manifold duration. 
%      - The manifold segment is constructed by propagating the specific 
%        trajectory in the halo orbit's stable manifold that corresponds to
%        the given tau value. The trajcetory depards the halo orbit either
%        in the interior or exterior direction as indicated. It is
%        popagated in the Earth Moon three body system backward in time for
%        the given duration.
% 3) Define XMI to be the final state of the manifold segment, This is the
%    state that a spacecraft would need to obtain in order to inject onto the
%    manifold segment, 
% 4) Construct DeltaXMI and the bridge segment. 
%      - Define DeltaXMI to be the tangential deltaV that may be applied to
%      XMI in order to construct the bridge segment. 
%      - When propagated further backward in time, the bridge segment will
%      encounter the propagated 185 km LEO orbit at the bridge's first
%      perigee point. The bridge segment is propagated in the Earth-Moon
%      three body system. 
% 5) Construct DeltaVLEO, the tangential DeltaV that may be applied to
%    transfer the spacecraft from its LEO orbit onto the bridge segment. 

startpoints = 41;
orbitPeriod =  t_half(end)*2;
clear t_halo X_halo
% Generate new matrix for only specified values in the orbit. 41 points
% evenly distributed around the orbit are picked. 
[t_halo,X_halo] = ode45(@CRTBPmodel, linspace(0,orbitPeriod,startpoints), X0, ode_options,mu);

% Get the monodromy matrix (the STM propagated     
% for an entire orbit. M = STM(t0+P,t0)), The matrix contains information
% aobtut every region that a spacecraft would pass through along that
% orbit. 
M = reshape(X_halo(end,7:end),6,[]); 

[eigVec eigVal] = eig(M); % eigenvalues and vectors of the monodromy matrix.
cnt = 0;
% Find the eigenvector of the stable eigenvector. (real and not maximum)
egenvarde = diag(eigVal);
for inda = 1:length(eigVal)
    tf = isreal(egenvarde(inda));
    if tf ==1
        cnt = cnt+1;
        values(cnt) = egenvarde(inda);
    end
end
[vv,place] = max(values);
[vv1,place1] = min(values);
vS = eigVec(:,place1);
vU = eigVec(:,place);

%Create starting points for the manifold. 
farg = get(gca,'colororder');
bromsrop = 0;
figure(4)
hold on
for indb = 1:startpoints
    % Make the pertubation vector.
	epsilonx =100/EM_dist; % Small pertubation, how small?
    epsilonv = epsilonx/norm(X_halo(indb,1:3));% Perturbation in velocity.
    epsilon = [epsilonx,epsilonx,epsilonx,epsilonv,epsilonv,epsilonv]';

    % The stable Manifold may be mapped by propagating the state XS backwards
    % in time where XS = X+-epsilon*vS.
    % Map to end
    % Use STM at correct time to map to the right time. 
    vStable = reshape(X_halo(indb,7:end),6,[])*vS;
    vStable = (epsilon.*vStable)/norm(vStable);
    
    XSInterior(indb,:) = X_halo(indb,1:6)'+vStable;% Add the pertubation. Plus or minus for exterior or interior manifold.
    XSExterior(indb,:) = X_halo(indb,1:6)'-vStable; %Exterior Stable
   
    % Propagate the manifold, for each of the starting points using the three
    % body system. 
    ode_options2 = odeset('Events',@FindLEO,'RelTol',1e-13,'AbsTol',1e-13);
    [t_manifoldI,X_manifoldI] = ode45(@CRTBPmodelNOSTM, [0 -Inf], XSInterior(indb,:), ode_options2,mu);
    [t_manifoldE,X_manifoldE] = ode45(@CRTBPmodelNOSTM, [0 -0.90*pi], XSExterior(indb,:), ode_options,mu);
    figure(4)
    plot3(X_manifoldI(:,1),X_manifoldI(:,2),X_manifoldI(:,3),'color',farg(4,:));hold on
    grid on
    figure(4)
    plot3(X_manifoldE(:,1),X_manifoldE(:,2),X_manifoldE(:,3),'color',farg(5,:));hold on
    pause(0.001)
    XMI = X_manifoldI(end,1:6);
    % Now We need to construct the deltaVMi and the bridge segment. 
    % Find the change in velocity that takes us to a 185km parking orbit at
    % Earth.

    ode_options2 = odeset('Events',@FindLEO,'RelTol',1e-13,'AbsTol',1e-13);
%     circle(-mu,0,parkingorbit_Earth);
    for indg = linspace(0,100,100)
        XMI1(1) = XMI(1);
        XMI1(2) = XMI(2);
        XMI1(3) = XMI(3);
        XMI1(4) =  XMI(4) + XMI(1)/norm([XMI(1:3)])*indg;
        XMI1(5) =  XMI(5) + XMI(2)/norm([XMI(1:3)])*indg;
        XMI1(6) =  XMI(6) + XMI(3)/norm([XMI(1:3)])*indg;
        [t_MI,X_MI] = ode45(@CRTBPmodelNOSTM, [0 -2*pi], XMI1, ode_options2,mu);

         pause(0.001)
         view(0, 90)
        for indc = 1:length(X_MI(:,1))
            if ( norm(X_MI(indc,1:3))<mu )&& (((X_MI(indc,2))<0.001) && (X_MI(indc,2)>-0.001)) %((X_MI(indc,1))>-(mu-0.001)||((X_MI(indc,1))<-(mu+0.001)))&&(((X_MI(indc,2))<0.001)||((X_MI(indc,2))>-0.001))
                bromsrop = 1;
                found = indc;
                break
            end
        end
        if bromsrop ==1;
            bromsrop =0;
            break
        end
    end
    plot3(X_MI(1:found,1),X_MI(1:found,2),X_MI(1:found,3),'m');hold on
end

grid on
figure(4)
plot3(0-mu,0,0,'ok', 'markerfacecolor', 'b', 'markersize', 22);hold on
plot3(1-mu,0,0,'ok', 'markerfacecolor', 'y', 'markersize',10);hold on
plot3(0.8369151324,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 8);hold on
plot3(1.1556821603,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 8);hold on
view(0, 90)
figure(10)

plot3(0-mu,0,0,'ok', 'markerfacecolor', 'b', 'markersize', 22);hold on
plot3(1-mu,0,0,'ok', 'markerfacecolor', 'y', 'markersize',10);hold on
plot3(0.8369151324,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 8);hold on
plot3(1.1556821603,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 8);hold on
view(0, 90)
plot3(X_manifoldI(:,1),X_manifoldI(:,2),X_manifoldI(:,3),'g','linewidth',1.5);hold on
plot3(X_MI(1:found,1),X_MI(1:found,2),X_MI(1:found,3),'m','linewidth',1.5);hold on
plot3((X_halo(:,1)),X_halo(:,2),X_halo(:,3),'b','linewidth',1.5);hold on
title('Transfer to Halo');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
grid on