function [ h ] = WeibulProcessor( lifetime, datamode )
%WeibulProcessor Generate hazard curve on basis of selected statistical data
%   Detailed explanation goes here
%% Weibull setup
datamode = string(datamode); %micro, 2000, comms, comp 

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

% labelloc = 0:365:365*lifetime;
% plot(t,R_micro,t,R_2000,t,R_comms,t,R_comp,'LineWidth',1.5);
% legend('Microsat','After 2000','comms','comp')
% legend('prob density','prob dist','Reliability')
% plot(t,f,t,F,t,R,'LineWidth',1.5)
% 
% 
% xlim([0 t(end)])
% xticks(labelloc)
% xticklabels({'0','1','2','3','4','5'})
% xlabel('Time [years]')
% ylabel('Probability [-]')


end
