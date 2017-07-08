%% input parameters
n = 6; % orbit planes
p = 6; % sats/orbit
lifetime = 5; %years
n_cycles = 10000;

%Request weibulcurve
[h,t] = WeibulProcessor(lifetime,'newer',0);
prob = 1-h;
T_replace = 10; % time to replace a broken satellite [days]
%t = 1:365*lifetime;

nsamples = t(end);

%Construct satellite IDs
orbits = 100:100:n*100;
sats = 1:1:p;
ID = repmat(orbits,p,1) + repmat(sats',1,n);
[rdim,cdim] = size(ID);

%Initialize
maxbroken = [];
histlist = [];
sysfailtime = zeros(t(end),n_cycles);
brokenID = [];
faillist = zeros(t(end),n_cycles);

for cycles=1:n_cycles
    satfail = 0;
    sysfail = 0;
    failures = zeros(size(ID));
    status = zeros(size(ID));
    rows = [];
    cols = [];
    failure_IDs = [];
    stats = [];
    broken = [];
    
    if cycles/n_cycles == 0.10||cycles/n_cycles == 0.20||cycles/n_cycles == 0.30||cycles/n_cycles == 0.40||cycles/n_cycles == 0.50||cycles/n_cycles == 0.60||cycles/n_cycles == 0.70||cycles/n_cycles == 0.80||cycles/n_cycles == 0.90
        disp(['%'])
    end

    for s=1:nsamples
        faillist(s,cycles) = satfail;
        replaced = [];
        sample = rand(size(ID));%exprnd(0.0950,size(ID)); % get random sample

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
                    notice_satfail = ['Sats ' num2str(failure_ID') ' have failed. '];
                    %disp(notice_satfail)

                    stats = [rows cols failure_IDs]';
                end
            end
        end

        %Loop through stats register
        for l=1:size(stats,2)
            e = stats(:,l);
            %Repair progress and spares check

            if status(e(1),e(2))>=T_replace %sat is replaced

                %notice_repair = ['Sat ' num2str(e(3)) ' is replaced.'];
                %disp(notice_repair)
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
                if abs(diff)==100||abs(diff)==100*(size(ID,1)-1)%||abs(diff)==1||abs(diff)==size(ID,2)-1
                    %notice_sysfail = ['Adjacent sat failure. ' num2str(failure_IDs') ' have failed simultaneously.'];
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

%% plot
flist = mean(faillist,2);
downfrac = mean(histlist);
timefrac = mean(sum(sysfailtime))/(365*lifetime);
areafrac = timefrac*2/(n*p);
notice_end = ['Adjacent failure%: ' num2str(downfrac*100) ', downtime%: ' num2str(timefrac*100) ', Total broken: ' num2str(flist(end))];
disp(notice_end)
figure;
yyaxis left

plot(t,mean(sysfailtime,2),'LineWidth',1.5)
ylabel('Adjacent plane failure probabilty','FontSize',12)
ylim([0 0.004])
yyaxis right
flist = mean(faillist,2);
plot(1:length(flist),flist,'LineWidth',1.5)
ylabel('Total broken satellites','FontSize',12)
labelloc = 0:365:365*lifetime;
xlim([0 t(end)])
xticks(labelloc)
xticklabels(0:lifetime)
xlabel('Time [years]','FontSize',12)

% dim = [.2 .3 0 0];
% str = 'Straight Line Plot from 1 to 10';
% annotation('textbox',dim,'String',notice_end,'FitBoxToText','on');