%inputs
p = 6; % sats/orbit
n = 1; % single orbit analysis
lifetime = 5; %years
n_cycles = 1000;

T_resend = 183; % days in 6 months of resend time through launch
T_replace = 10; % time to replace a broken satellite from spares

%t = 1:365*lifetime;
[h,t] = WeibulProcessor(lifetime,'new',0);
prob = 1-h;

%launch settings
n_spare_start = 4; %spares at initial condition
launch_thresh = 1; %spares amount when new launch is ordered
n_restock = 2; %amount of spares added to a plane at launch arrival

%Construct satellite IDs
orbits = 100:100:n*100;
sats = 1:1:p;
ID = repmat(orbits,p,1) + repmat(sats',1,n);
[rdim,cdim] = size(ID);

iterator = 1:4;

dat1 = zeros(length(iterator),1);
dat2 = zeros(length(iterator),1);
dat3 = zeros(length(iterator),1);

lable = 'n_restock'; %EDIT THIS

for itvl=iterator
    n_restock = itvl; %EDIT THIS
    disp(itvl);
    %Initialise
    maxbroken = [];
    histlist = [];
    sysfailtime = zeros(t(end),n_cycles);
    brokenID = [];
    launchlist = zeros(t(end),n_cycles);
    faillist = zeros(t(end),n_cycles);
    sparelist = zeros(t(end),n_cycles);
    nsamples = t(end);
    
    for cycles=1:n_cycles
        %initialise
        n_l = 0;
        satfail = 0;
        sysfail = 0;
        failures = zeros(size(ID));
        status = zeros(size(ID));
        rows = [];
        cols = [];
        failure_IDs = [];
        stats = [];
        broken = [];
        launchnew = 0;
        n_spare = n_spare_start;
        Tsend = 0;

        for s=1:nsamples
            sparelist(s,cycles) = n_spare;
            launchlist(s,cycles) = n_l;
            faillist(s,cycles) = satfail;
            replaced = [];
            sample = rand(size(ID));%exprnd(0.0950,size(ID)); % get random sample
            if launchnew == 1 %if a launch is being conducted, update launch arrival timer
                Tsend = Tsend + 1;
                if Tsend >= T_resend %if launch arrived, update spares
                    n_spare = n_spare + n_restock;
                    Tsend = 0;
                    launchnew = 0;
                    n_l = n_l + 1;
                    %disp('LAUNCH ARRIVED')
                end

            end

            if any(any(sample>prob(s))) % if initial failure occurs
                failmat = ID.*(sample>prob(s));
                failure_ID = (failmat(failmat~=0));

                [row,col] = find(failmat);

                %Looping through just failed satellites
                for f=1:length(failure_ID)

                    %Check if satellite already failed
                    if any(failure_ID(f)==failure_IDs)
                        %Do nothing, satellite has already failed
                    else %Add failed satellite to stats register
                        satfail = satfail + 1;
                        rows = [rows;row];
                        cols = [cols;col];
                        failure_IDs = [failure_IDs;failure_ID];
                        notice_satfail = ['Sats ' num2str(failure_ID') ' have failed.'];
                        %disp(notice_satfail)

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
                        %disp(notice_repair)
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
                        notice_sysfail = ['Adjacent sat failure. ' num2str(failure_IDs') ' have failed simultaneously. Spares left: ' num2str(n_spare)];
                        sysfail = 1;
                        sysfailtime(s,cycles) = 1;
                        bID = failure_IDs;

                    end
                end
            end  

        end
        if sysfail == 1
            %disp(notice_sysfail)
            brokenID = [brokenID bID'];

        end
        maxbroken = [maxbroken max(broken)];
        histlist = [histlist sysfail];
    end
    
%     slist = mean(sparelist,2);
%     llist = mean(launchlist,2);
%     flist = mean(faillist,2);

% 
%     endfail = faillist(1825,:);
     endlaunch = launchlist(1825,:);
     minspare = min(sparelist);
     sftlist = mean(sysfailtime,2);

     dat1(itvl-iterator(1)+1) = sum(minspare==0)/n_cycles;
     dat2(itvl-iterator(1)+1) = mean(endlaunch);
     dat3(itvl-iterator(1)+1) = sum(sftlist)/t(end);
end

figure;
xlabel(lable)
yyaxis left
plot(iterator,dat1);
ylabel('mean days of 0 spares present')


yyaxis right
plot(iterator,dat2);
ylabel('mean launches after 5 years')
figure;
plot(iterator,dat3);
xlabel(lable)
ylabel('lifetime system failure fraction')

