%input parameters
n = 4; % orbit planes
p = 8; % sats/orbit
%ASSUME 1 SAT FAILS PER YEAR
nsf = 1/(n*p); %#sats failing per year
prob = 1-nsf/(365*24); %1- #sats failing per hour = prob
T_replace = 365*24/2; %hours in 6 months of replace time
nsamples = 131400; %365*24*15

%Construct satellite IDs
orbits = 100:100:n*100;
sats = 1:1:p;
ID = repmat(orbits,p,1) + repmat(sats',1,n);
[rdim,cdim] = size(ID);

%Initialize
failures = zeros(size(ID));
status = zeros(size(ID));
rows = [];
cols = [];
failure_IDs = [];
stats = [];
broken = [];
sysfail = 0;
maxbroken = [];
count = 0;
histlist = [];
failtime = 0;
brokenID = [];
for cycles=1:50
    sysfail = 0;
    count = 0;
    while count < 1&&sysfail==0
        disp(count)
        failures = zeros(size(ID));
        status = zeros(size(ID));
        rows = [];
        cols = [];
        failure_IDs = [];
        stats = [];
        broken = [];

        for s=1:nsamples
            replaced = [];
            sample = rand(size(ID)); % get random sample

            if any(any(sample>prob)) % if initial failure occurs
                failmat = ID.*(sample>prob);
                failure_ID = (failmat(failmat~=0));

                [row,col] = find(failmat);

                %Looping through just failed satellites
                for f=1:length(failure_ID)

                    %Check if satellite already failed
                    if any(failure_ID(f)==failure_IDs)
                        %Do nothing, satellite has already failed
                    else %Add failed satellite to stats register
                        rows = [rows;row];
                        cols = [cols;col];
                        failure_IDs = [failure_IDs;failure_ID];

                        stats = [rows cols failure_IDs]';
                    end
                end
            end

            %Loop through stats register
            for l=1:size(stats,2)
                e = stats(:,l);
                %Repair progress
                if status(e(1),e(2))>=T_replace %sat is replaced
                    notice_repair = ['Sat ' num2str(e(3)) ' is replaced.'];
                    %disp(notice_repair)
                    %stats(:,l)=[];
                    replaced = [replaced l];
                    status(e(1),e(2)) = 0;
                else
                    status(e(1),e(2)) = status(e(1),e(2)) + 1;
                end
            end
            %Remove repaired sats from the stat register
            %stats(:,replaced)=[];
            rows(replaced) = [];
            cols(replaced) = [];
            failure_IDs(replaced) = [];
            stats(:,replaced) = [];
            broken = [broken size(stats,2)];

            %Analyse adjacent failure possibilities
            if broken(length(broken))>1 %if more than 2 are not operational simultaneously
                IDdiff = [];
                %generate ID difference list
                for br=1:(length(failure_IDs)-1)
                    for brr=(br+1):(length(failure_IDs))
                        IDdiff = [IDdiff failure_IDs(br)-failure_IDs(brr)];
                    end
                end
                

                for diff=IDdiff
                    if abs(diff)==100||abs(diff)==1||abs(diff)==11||abs(diff)==500
                        notice_sysfail = ['Adjacent sat failure. ' num2str(failure_IDs') ' have failed simultaneously.'];
                        sysfail = 1;
                        failtime = failtime + 1;
                        bID = failure_IDs;
                        
                    end
                end
            end     
        end
        if sysfail == 1
            disp(notice_sysfail)
            brokenID = [brokenID bID'];
            
        end
        maxbroken = [maxbroken max(broken)];
        count = count + 1;

    end
    histlist = [histlist count*sysfail];
end
notice_end = ['Total non-100% time [h]: ' num2str(failtime) '. Max sats broken simultaneously: ' num2str(max(maxbroken))];
disp(notice_end)
hist(histlist)
