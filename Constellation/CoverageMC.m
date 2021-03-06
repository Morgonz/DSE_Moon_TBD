%% Weibull setup
datamode = 'comms'; %micro, 2000,comms,comp 
lifetime = 5;% years
t = 1:365*lifetime; %Time in days
%Weibull for microsats:
beta =  0.2928; 
eta = 10065; % days
R_micro = exp(-(t/eta).^beta);
%Weibull for small sats after 2000:
beta =  0.3256; 
eta = 2180; % days 
R_2000 = exp(-(t/eta).^beta);
%Weibull for communication missions:
beta =  0.4023; 
eta = 6524; % days
R_comms = exp(-(t/eta).^beta);
%Weibull for commercially produced:
beta =  0.3318; 
eta = 1629; % days
R_comp = exp(-(t/eta).^beta);

if strcmp('micro',datamode)
    %Weibull for microsats:
    beta =  0.2928; 
    eta = 10065; % days
    disp('datamode = micro')
elseif strcmp('2000',datamode)
    %Weibull for small sats after 2000:
    beta =  0.3256; 
    eta = 2180; % days 
    disp('datamode = 2000')
elseif strcmp('comp',datamode)
    %Weibull for communication missions:
    beta =  0.4023; 
    eta = 6524; % days
    disp('datamode = comms')
elseif strcmp('comms',datamode)
    %Weibull for commercially produced:
    beta =  0.3318; 
    eta = 1629; % days
    disp('datamode = comp')
else
    disp('datamode not recognised')
end
R = exp(-(t/eta).^beta); %Set probability distribution with chosen dataset
F = 1-R;
f = beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta);
h = f./R;
labelloc = 0:365:365*lifetime;

plot(t,R_micro,t,R_2000,t,R_comms,t,R_comp,'LineWidth',1.5);
legend('Microsat','After 2000','comms','comp')
% legend('prob density','prob dist','Reliability')
% plot(t,f,t,F,t,R,'LineWidth',1.5)
% 
% 
% xlim([0 t(end)])
% xticks(labelloc)  
% xticklabels({'0','1','2','3','4','5'})
% xlabel('Time [years]')
% ylabel('Probability [-]')







%% input parameters
n = 6; % orbit planes
p = 6; % sats/orbit

T_resend = 183; %days in 6 months of resend time
T_replace = 10; %one month

launchnew = 0;
n_spare_start = 24;
launch_thresh = 8;
n_restock = 12;
nsamples = t(end);
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
launchlist = [];
n_cycles = 1;
faillist = [];
sparelist = [];
n_l = 0;
satfail = 0;
for cycles=1:n_cycles
    sysfail = 0;
    count = 0;

    while count < 1&&sysfail==0
        st = [cycles count];
        %disp(st)
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
            prob = 1-h(s);
            sparelist = [sparelist n_spare];
            launchlist = [launchlist n_l];
            faillist = [faillist satfail];
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
            prob

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
                        satfail = satfail + 1;
                        rows = [rows;row];
                        cols = [cols;col];
                        failure_IDs = [failure_IDs;failure_ID];
                        notice_satfail = ['Sats ' num2str(failure_ID') ' have failed. '];
                        disp(notice_satfail)

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
notice_end = ['Total non-100% time%: ' num2str(failtime/(365*15)*100*sum(histlist)/n_cycles) '. Area% not covered: ' num2str(failtime/(365*15*n*p)*2*100*sum(histlist)/n_cycles) '. Max sats broken simultaneously: ' num2str(max(maxbroken))];
disp(notice_end)
figure;
hist(histlist)
figure(2);
yyaxis left
plot(1:length(sparelist),sparelist,'LineWidth',1.5)
title('Spares status','FontSize',16)
xlabel('Time [days]','FontSize',12)
ylabel('Spares remaining','FontSize',12)

yyaxis right
plot(1:length(faillist),faillist,'LineWidth',1.5)
ylabel('Total broken satellites','FontSize',12)
xlim([0 length(faillist)])