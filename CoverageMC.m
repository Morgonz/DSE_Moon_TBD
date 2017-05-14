%Initialise
n = 6; % orbit planes
p = 12; % sats/orbit
prob = 0.9; %Probability of one satellite failing per HOUR
T_replace = 365*24/2; %hours in 6 months of replace time
nsamples = 100;

%Construct satellite IDs
orbits = 100:100:n*100;
sats = 1:1:p;
ID = repmat(orbits,p,1) + repmat(sats',1,n);

%ID = [101:112; 201:212; 301:312; 401:412; 501:512; 601:612];
failures = zeros(size(ID));
for i=1:nsamples
    sample = rand(size(ID));
    if any(any(sample>prob))
        failure = ID.*(sample>prob);
        failure_ID = failure(failure~=0);
        for j=1:T_replace
            %time to repair
        end
    end
    
    %else: no failure...
    
end
out = failures/nsamples



% failures = zeros(size(ID));
% for i=1:nsamples
%     failures = failures + (rand(size(ID))>prob);
%     
% end
% out = failures/nsamples