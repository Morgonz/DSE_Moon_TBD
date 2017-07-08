function [ h,t ] = WeibulProcessor( lifetime, datamode, plotall )
%WeibulProcessor Generate hazard curve on basis of selected statistical data
%   Lifetime in years, datamode inputs are: micro, 2000, comms, comp, new

%% Weibul data selection
datamode = string(datamode); %micro, 2000, comms, comp, new

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
%New weibull distribution
alpha = 0.07102;
eta1 = 368400*365;
eta2 = 6.492*365;
beta1 = 0.05344;
beta2 = 6.859;
R_new = alpha*exp(-(t/eta1).^beta1)+(1-alpha)*exp(-(t/eta2).^beta2);
%Newer weibull distribution
alpha = 0.9607;
eta1 = 10e7*365;
eta2 = 7.3*365;
beta1 = 0.2101;
beta2 = 2.754;
R_3 = alpha*exp(-(t/eta1).^beta1)+(1-alpha)*exp(-(t/eta2).^beta2);
f_3 = alpha*beta1/eta1.*(t/eta1).^(beta1-1).*exp(-(t/eta1).^beta1)+(1-alpha)*beta2/eta2.*(t/eta2).^(beta2-1).*exp(-(t/eta2).^beta2);
h_3 = f_3./R_3;

if strcmp('micro',datamode)
    %Weibull for microsats:
    beta =  0.2928; 
    eta = 10065; % days
    disp('datamode = micro')
    R = exp(-(t/eta).^beta); %Set probability distribution with chosen dataset
    f = beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta);
    h = f./R;
elseif strcmp('2000',datamode)
    %Weibull for small sats after 2000:
    beta =  0.3256; 
    eta = 2180; % days 
    disp('datamode = 2000')
    R = exp(-(t/eta).^beta); %Set probability distribution with chosen dataset
    f = beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta);
    h = f./R;
elseif strcmp('comp',datamode)
    %Weibull for communication missions:
    beta =  0.4023; 
    eta = 6524; % days
    disp('datamode = comms')
    R = exp(-(t/eta).^beta); %Set probability distribution with chosen dataset
    f = beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta);
    h = f./R;
elseif strcmp('comms',datamode)
    %Weibull for commercially produced:
    beta =  0.3318; 
    eta = 1629; % days
    disp('datamode = comp')
    R = exp(-(t/eta).^beta); %Set probability distribution with chosen dataset
    f = beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta);
    h = f./R;
elseif strcmp('manual',datamode)
    %Manually adapted Weibul:
    beta =  0.1; 
    eta = 1629; % days
    disp('datamode = comp')
    R = exp(-(t/eta).^beta); %Set probability distribution with chosen dataset
    f = beta/eta.*(t/eta).^(beta-1).*exp(-(t/eta).^beta);
    h = f./R;
elseif strcmp('new',datamode)
    alpha = 0.07102;
    eta1 = 368400*365;
    eta2 = 6.492*365;
    beta1 = 0.05344;
    beta2 = 6.859;
    R = alpha*exp(-(t/eta1).^beta1)+(1-alpha)*exp(-(t/eta2).^beta2);
    R = horzcat(R(1:1387),0.947*ones(1,t(end)-1387));
    f = alpha*beta1/eta1.*(t/eta1).^(beta1-1).*exp(-(t/eta1).^beta1)+(1-alpha)*beta2/eta2.*(t/eta2).^(beta2-1).*exp(-(t/eta2).^beta2);
    h = f(1:1387)./R(1:1387);
    h_2end = h(end);
    h = horzcat(h,h_2end*ones(1,t(end)-1387));
%     alpha = 0.07102;
%     eta1 = 368400*365;
%     eta2 = 6.492*365;
%     beta1 = 0.05344;
%     beta2 = 6.859;
%     R = alpha*exp(-(t/eta1).^beta1)+(1-alpha)*exp(-(t/eta2).^beta2);
%     f = alpha*beta1/eta1.*(t/eta1).^(beta1-1).*exp(-(t/eta1).^beta1)+(1-alpha)*beta2/eta2.*(t/eta2).^(beta2-1).*exp(-(t/eta2).^beta2);
%     h = f./R;
elseif strcmp('newer',datamode)
    alpha = 0.9607;
    eta1 = 1e7*365;
    eta2 = 7.3*365;
    beta1 = 0.2101;
    beta2 = 2.754;
    R = alpha*exp(-(t/eta1).^beta1)+(1-alpha)*exp(-(t/eta2).^beta2);
    f = alpha*beta1/eta1.*(t/eta1).^(beta1-1).*exp(-(t/eta1).^beta1)+(1-alpha)*beta2/eta2.*(t/eta2).^(beta2-1).*exp(-(t/eta2).^beta2);
    h = f./R;
else
    error('datamode not recognised. Please use "micro", "2000", "comms", "comp" or "new" as datamode input.')
    
end

%% Generate Weibul curve and determine hazard rate


if plotall
    figure;
    plot(t,R_micro,t,R_2000,t,R_comms,t,R_comp,t,R_new,t,R_3,'LineWidth',1.5);
    legend('Microsat','After 2000','comms','comp','new','newer')
    labelloc = 0:365:365*lifetime;
    axis([0 t(end) 0 1])
    xticks(labelloc)
    xticklabels(0:lifetime)
    xlabel('Time [years]')
    ylabel('Probability [-]')
    title('Probability distribution overview')
    
    figure;    
    yyaxis left
    plot(t,R,'LineWidth',2);
%     plot(t,f,'LineWidth',2);hold off;
    ylim([R(end) 1])
    ylabel('Probability [-]','FontSize',12)
    
    yyaxis right
    plot(t,1-h,'LineWidth',2)
    ylim([1-h(1) 1])
    ylabel('Failure rate [-]','FontSize',12)

    labelloc = 0:365:365*lifetime;
    xlim([0 t(end)])
    xticks(labelloc)
    xticklabels(0:lifetime)
    xlabel('Time [years]','FontSize',12)
%     title('Selected mode details')
    lgnd = legend('2-Weibull mixture reliability curve','Hazard curve','Location','southwest');
    set(lgnd,'FontSize',10); 
end
% Verification data
% points_h = [1 2 5 20 100 500 1800];
% result_h = [0.9802963267 0.9879313474 0.9936870026 0.9976315796 0.9992411653 0.9997568916 0.9999017306];
% points_R = [1 25 200 500 1500];
% result_R = [0.934920326 0.8413897857 0.727980529 0.6602226179 0.5639942852];
end