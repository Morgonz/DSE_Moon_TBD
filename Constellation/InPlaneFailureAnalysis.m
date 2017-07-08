%inputs
p = 6; % sats/orbit
n = 1; % single orbit analysis
lifetime = 15; %years
n_cycles = 10000;

%t = 1:365*lifetime;
[h,t] = WeibulProcessor(lifetime,'newer',0);
prob = 1-h;

%launch settings
T_resend = 182; % days in 6 months of resend time through launch
T_replace = 10; % time to replace a broken satellite from spares

n_spare_start = 2; %spares at initial condition
launch_thresh = 0; %spares amount when new launch is ordered
n_restock = 2; %amount of spares added to a plane at launch arrival

%Construct satellite IDs
orbits = 100:100:n*100;
sats = 1:1:p;
ID = repmat(orbits,p,1) + repmat(sats',1,n);
[rdim,cdim] = size(ID);

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
    
    
    if cycles/n_cycles == 0.10||cycles/n_cycles == 0.20||cycles/n_cycles == 0.30||cycles/n_cycles == 0.40||cycles/n_cycles == 0.50||cycles/n_cycles == 0.60||cycles/n_cycles == 0.70||cycles/n_cycles == 0.80||cycles/n_cycles == 0.90
        disp(['%'])
    end

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
            else 
                launchnew = 1;
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
%% Plotting
%downfrac = mean(histlist);
%timefrac = sysfailtime/(365*lifetime*n_cycles);
%areafrac = timefrac*2/(n*p);
%notice_spares = ['initial amount: ' num2str(n_spare_start) ' l_thresh: ' num2str(launch_thresh) ' n_restock: ' num2str(n_restock)];
%disp(notice_spares)
%notice_end = ['Total cycle sysfail%: ' num2str(downfrac*100) ' Total non-100% time%: ' num2str(timefrac*100) '. Area% not covered: ' num2str(areafrac*100)];
%disp(notice_end)


% x = sparelist;
% xu = mean(x,2) + norminv(0.95) * (std(sparelist,0,2) / (size(sparelist,2))^2) ;
% xl = mean(x,2) - norminv(0.95) * (std(sparelist,0,2) / (size(sparelist,2))^2) ;

figure; %Data means with confidence intervals TBD
yyaxis right
slist = mean(sparelist,2);
plot(1:length(slist),slist,'LineWidth',1); hold on;
llist = mean(launchlist,2);
plot(1:length(llist),llist,'LineWidth',1); 
flist = mean(faillist,2);
plot(1:length(flist),flist,'LineWidth',1.5); hold off;
% title('Spares status per orbit: Averages','FontSize',16)
xlabel('Time [days]','FontSize',12)
ylabel('Number','FontSize',12)
ylim([0 1.25*max([max(slist),max(flist),max(llist)])])

yyaxis left
sftlist = mean(sysfailtime,2);
plot(1:length(sftlist),sftlist,'LineWidth',1); 
ylabel('Plane system failure probability [-]','FontSize',12)
ylim([0 max(sftlist)*1.1+0.00005])

xlabelloc = 0:365:365*lifetime;
xlim([0 t(end)])
xticks(xlabelloc)
xticklabels(0:lifetime)
xlabel('Time [years]','FontSize',12)
legend('Adjacent failure probability','Average spares in plane','Average launches arriving','Average broken satellites','Location','northeast')
set(gca,'FontSize',12)
figure; %Case studies
% casetitle = suptitle('Statistical analysis');
% set(casetitle,'FontSize',18,'FontWeight','normal')

subplot(1,3,1) %Satellite failures plot
endfail = faillist(1825,:);
histogram(endfail,'EdgeAlpha',0.4,'FaceColor','r','Normalization','probability'); 
xlabel('Failed sats at 5 years','FontSize',12)
ylabel('Probability of occurrence [-]','FontSize',14)
xlim([-0.5 max(endfail)+0.5])
xticks(0:max(endfail))
xticklabels(0:max(endfail))
set(gca,'FontSize',12)
subplot(1,3,2) % Launch amount plot
endlaunch = launchlist(1825,:);
histogram(endlaunch,'EdgeAlpha',0.4,'FaceColor','b','Normalization','probability');
xlabel('Launches at 5 years','FontSize',12)
xlim([-0.5 max(endlaunch)+0.5])
xticks(0:max(endlaunch))
xticklabels(0:max(endlaunch))
set(gca,'FontSize',12)
subplot(1,3,3) % Minimum spares present plot
minspare = min(sparelist);
histogram(minspare,'EdgeAlpha',0.4,'FaceColor','g','Normalization','probability');
xlabel('Minimum spares in plane','FontSize',12)
xlim([-0.5 max(minspare)+0.5])
xticks(0:max(minspare))
xticklabels(0:max(minspare))
set(gca,'FontSize',12)
notice_out = ['Spares in orbit: ' num2str(slist(end)) '. Launches: ' num2str(llist(end)) '. Sat failures: ' num2str(flist(end)) '. Prob of a failure in lifetime: ' num2str(mean(mean(sysfailtime)>0)*100) '. Downtime frac:' num2str(mean(mean(sysfailtime))*100) '.'];
disp(notice_out)