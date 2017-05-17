%% Visibility checker


%% Create Moon Mesh
Rm  = 1738100; % [m]

[X,Y,Z] = sphere(1000);
X = reshape(X,length(X)*length(X),1).*Rm;
Y = reshape(Y,length(Y)*length(Y),1).*Rm;
Z = reshape(Z,length(Z)*length(Z),1).*Rm;

Moon_mesh = [X,Y,Z]; % Npoints x 3

%% Import Satellite locations % Nsats x 3 
sat_loc =    [600e3     0        0; 
                0      600e3     0;
                0,      0      600e3];
sat_loc = (sat_loc~=0)*Rm+sat_loc;

%% Check Visibility
max_dist = 1200e3; %m
vis_list = zeros(1,length(Moon_mesh)); %empty list to store the jdx of the visible ground locations
%loop through the ground mesh
for idx=1:length(sat_loc)
    for jdx=1:length(Moon_mesh)
        u_L = norm(Moon_mesh(jdx,:)-sat_loc(idx,:));
        if u_L<=max_dist
            vis_list(jdx) = 1;
        end   
    end  
end

figure
scatter3(X,Y,Z,1,vis_list)

    