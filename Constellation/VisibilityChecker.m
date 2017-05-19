%% Visibility checker
%clear all;
close all;
%clc;

%% Create Moon Mesh
Rm  = 1738100+700e3; % [m]

[X,Y,Z] = sphere(200);
X = reshape(X,length(X)*length(X),1).*Rm;
Y = reshape(Y,length(Y)*length(Y),1).*Rm;
Z = reshape(Z,length(Z)*length(Z),1).*Rm;

Moon_mesh = [X,Y,Z]; % Npoints x 3

%% Import Satellite locations % Nsats x 3 
% Network
%sat_loc =    [A(1)     B(1)        C(1); 
            %  D(1)     E(1)        F(1)];
%sat_loc = (sat_loc~=0)*Rm+sat_loc;

%sat_loc = walkerdeltagenerator(60,2,3,400e3);
sat_loc = loc_sat(:,1:3)

relay = 1;
% Relay
rel_loc =    [5000e3     0         0; 
                0      5000e3     0;
                0,      0      5000e3;
                0,      0      5000e3];
rel_loc = (rel_loc~=0)*Rm+rel_loc;


%% Set satellite and rover parameters - validate using existing sats

ele_ang = deg2rad(5); % [rad]

% Network
S_rec   = deg2rad(20); % [rad] half beam angle when satellite receives data from rover
max_dist_G_rec = 12000e3; % (norm(sat_loc(1,:))-Rm)/cos(S_rec); % [m] Max distance at S_rec to close budget

S_tra   = deg2rad(20); % [rad] half beam angle when satellite sends data to rover
max_dist_S_rec = 12000e3; % (norm(sat_loc(1,:))-Rm)/cos(S_tra); % [m] Max distance at S_rec to close budget

% Relay
R_rec   = deg2rad(25); % [rad] half beam angle when relay receives data from network
max_dist_R_rec = 10000e3; % [m] Max distance at R_rec to close budget

R_tra   = deg2rad(25); % [rad] half beam angle when relay sends data to network
max_dist_R_tra = 10000e3; % [m] Max distance at R_tra to close budget



%% Check Visibility
SG_vis_list = zeros(1,length(Moon_mesh)); %empty list to store the jdx of the visible ground locations
GS_vis_list = zeros(1,length(Moon_mesh));
SR_vis_list = zeros(1,length(sat_loc));
RS_vis_list = zeros(1,length(sat_loc));

%angle = atan2(norm(cross(a,b)),dot(a,b));

% loop through all the satellites
for idx=1:length(sat_loc)
    u_S = sat_loc(idx,:); % define the position vector of the satellite
    if norm(u_S)>=(Rm+100e3) % check wether satellite is lower than 100 km
        for jdx=1:length(Moon_mesh)
            u_G = Moon_mesh(jdx,:); %define pos vector for rover pos
            u_L = u_G-u_S;    % define vector going from S to R
            if atan2(norm(cross(-u_G,-u_L)),dot(-u_G,-u_L))>=(pi/2 + ele_ang)% check wether signal is above set elevation
                % rover -> networksat
                len_u_L = norm(u_L);
                if len_u_L<=max_dist_S_rec
                    if atan2(norm(cross(-u_S,u_L)),dot(-u_S,u_L))<=S_rec
                        %disp('hit1')
                        GS_vis_list(jdx) = GS_vis_list(jdx) + 1;    
                    end 
                end
                % networksat -> rover
                if len_u_L<=max_dist_G_rec
                    if atan2(norm(cross(-u_S,u_L)),dot(-u_S,u_L))<=S_tra
                        SG_vis_list(jdx) = SG_vis_list(jdx) + 1;    
                    end 
                end
            end   
        end    
        if relay==0
            for kdx=1:length(rel_loc)
            u_R = rel_loc(kdx,:);
            u_P = u_S-u_R;
            if norm(u_R)>=(norm(u_S)+100e3)
                len_u_P = norm(u_P);
                % networksat -> relaysat
                if len_u_P<=max_dist_R_rec
                    if atan2(norm(cross(-u_R,u_P)),dot(-u_R,u_P))<=R_rec
                        SR_vis_list(idx) =+ 1;    
                    end 
                end
                % relaysat -> networksat
                if len_u_P<=max_dist_R_tra
                    if atan2(norm(cross(-u_R,u_P)),dot(-u_R,u_P))<=R_tra
                        RS_vis_list(idx) =+ 1;    
                    end 
                end
            else 
            disp('error : Relay satellite altitude is lower than network alt + 100 km')
            Rsatfaulty = sat_loc(idx,:)
            end 
        end  
        end
    else
        disp('error : Network satellite altitude is lower than 100 km')
        Nsatfaulty = sat_loc(idx,:)
    end 
    %idx
end

zeroes = sum(GS_vis_list(:)==0)/length(GS_vis_list)*100;
ones = sum(GS_vis_list(:)==1)/length(GS_vis_list)*100;
rest = sum(GS_vis_list(:)>1)/length(GS_vis_list)*100;
notice_fin = [num2str(zeroes) ' ' num2str(ones) ' ' num2str(rest) ' R_sat: ' num2str(norm(sat_loc(1,:))-Rm) ' max dist: ' num2str(max_dist_S_rec) ' angle: ' num2str(rad2deg(S_rec))];
disp(notice_fin)

   

%% Plotting
% figure
% scatter3(X,Y,Z,1,SG_vis_list)
% hold on
% scatter3(sat_loc(:,1),sat_loc(:,2),sat_loc(:,3),10, SR_vis_list)
% % if relay==0
% %     scatter3(rel_loc(:,1),rel_loc(:,2),rel_loc(:,3),20,[0,0,0,0])
% % end
% colorbar
% axis equal

figure
scatter3(X,Y,Z,1,GS_vis_list)
hold on
scatter3(sat_loc(:,1),sat_loc(:,2),sat_loc(:,3),10,RS_vis_list)
if relay==0
    scatter3(rel_loc(:,1),rel_loc(:,2),rel_loc(:,3),20,[0,0,0,0])
end 
colorbar
axis vis3d
axis equal

% [lat,lon,z] = ecef2lla(X,Y,Z);
% lat = rad2deg(lat);
% lon = rad2deg(lon)-180;
% 
% figure
% scatter(lon,lat,1,GS_vis_list)
% 
% figure
% scatter(lon,lat,1,SG_vis_list)
% 


    