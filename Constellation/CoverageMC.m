%input parameters
n = 6; % orbit planes
p = 12; % sats/orbit
%ASSUME 1 SAT FAILS PER YEAR
nsf = 1/(n*p); %#sats failing per year
prob = 1-nsf/(365*24); %1- #sats failing per hour = prob

T_resend = 365*24/2; %hours in 6 months of resend time
T_replace = 24*365/12; %one month
launchnew = 0;
n_spare = 6;
launch_thresh = 3;
n_launch = 6;
nsamples = 131400; %365*24*15
Tsend = 0;

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

for cycles=1:10
    sysfail = 0;
    count = 0;
    while count < 1&&sysfail==0
        st = [cycles count];
        disp(st)
        failures = zeros(size(ID));
        status = zeros(size(ID));
        rows = [];
        cols = [];
        failure_IDs = [];
        stats = [];
        broken = [];
        launchnew = 0;
        n_spare = 6;
        Tsend = 0;

        for s=1:nsamples
            replaced = [];
            sample = rand(size(ID)); % get random sample
            
            if launchnew == 1 %if a launch is being conducted, update launch arrival timer
                Tsend = Tsend + 1;
                if Tsend >= T_resend %if launch arrived, update spares
                    n_spare = n_spare + n_launch;
                    Tsend = 0;
                    launchnew = 0;
                    %disp('LAUNCH ARRIVED')
                end
                    
            end

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
                %Repair progress and spares check
                if n_spare>0
                    if status(e(1),e(2))>=T_replace %sat is replaced
                        n_spare = n_spare - 1; %reduced spare count
                        if n_spare <= launch_thresh %if spares reaches launch lvl
                            launchnew = 1;

                            %disp('INITIATE LAUNCH')
                        end
                        notice_repair = ['Sat ' num2str(e(3)) ' is replaced. ' num2str(n_spare) ' spares left.'];
                        disp(notice_repair)
                        replaced = [replaced l];
                        status(e(1),e(2)) = 0;
                    else
                        status(e(1),e(2)) = status(e(1),e(2)) + 1;
                    end
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
notice_end = ['Total non-100% fraction: ' num2str(failtime/(24*365*15)) '. Max sats broken simultaneously: ' num2str(max(maxbroken))];
disp(notice_end)
hist(histlist)
