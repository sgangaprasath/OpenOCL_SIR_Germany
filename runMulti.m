%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Optimal policies for mitigating pandemic costs %%%%%%%%%%
% Credits: S Ganga Prasath, Mattia Serra, Salem Mosleh, Vidya Raju %
% Contact: gangaprasath@seas.harvard.edu
% Reference: arXiv:XXXX.XXXX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OpenOCL ref: Koenemann, Jonas, et al. "OpenOCL–Open Optimal
% Control Library." (2017).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all; clear pathall;
clear classes; clear java; clear functions;

addpath('./OpenOCL/OpenOCL-master/'); % Path of OpenOCL package

global Beta % Beta parameter in SIR model
global Gamma % Recovery rate in SIR model
global Tmax % Total time of simulation
global Np % Total population in millions
global Icus % No of ICUs for total population
global SInit % Initial susceptible population
global IInit % Initial infected population
global conf % Structure with country specific parameters
global N1 % Young population in millions
global N2 % Old population in millions
global uM % Upper bound for value of controller
global Py % Probability of young infected requiring ICU
global Po % Probability of old infected requiring ICU
global IFRy % IFR for young people 
global IFRo   % IFR for old people
global AlpL % Importance of Life cost
global AlpE % Importance of Economic cost
global AlpS % Importance of Social cost
global padding % Safety factor on total number of ICUs

AlpL = 1;
AlpE = 1;
AlpS = 1;

padding = 0.8;

Py = 0.0076;
Po = 0.031;

IFRy = 0.001;
IFRo = 0.02;

%% Germany parameters 
%------------------------------------------------------------------------------------------
Icus = 0.00034*padding;
Tmax = 365;
Beta = 0.036;
Gamma = 0.16;
N1 = 57;
N2 = 23;

IInit = [134e-6; 40e-6];
SInit = [N1; N2] - IInit;
Np = N1 + N2;

conf = struct;
conf.Acost = [Py/Icus 0; 0 Po/Icus]; % Normalized A matrix in economic cost
conf.Bcost = [1 0; 0 1]; % Weighing B matrix in economic cost 

conf.N1 = N1;
conf.N2 = N2;
Nin = diag([1/N1,1/N2]);
c011 = 8.5; c012 = 2.3; c021 = 5.6; c022 = 1.5;
conf.co = [c011 c012; c021 c022]*Nin; % Open-loop Control matrix
%------------------------------------------------------------------------------------------
%%
uM = .85;
Correc = 1/uM;
cc11 = c011; cc12 = Correc*c012; 
cc21 = c021;
cc22 = Correc*c022;
conf.cc = [cc11 cc12; cc21 cc22]*Nin; % Contribution to Closed-loop Control matrix

[solution,times,problem] = sirMultidim(); % Function to evolve SIR Optimal control dynamics 

greenish = [0, 255, 0]./255;
reddish = [255 0 0]./255;
blueish = [0 0 255]./255;

tFin = times.states.value;
sFin = [solution.states.S1.value./Np; solution.states.S2.value./Np;];
iFin = [solution.states.I1.value./Np; solution.states.I2.value./Np;];
rFin = [solution.states.R1.value./Np; solution.states.R2.value./Np;];
uFin = solution.controls.F.value;
iTotReal = sum(iFin);
iFiny = iFin(1,:);
iFino = iFin(2,:);
CostriFin = (Py*iFiny + Po*iFino)/Icus;
iTot = solution.states.It.value./Np;

% % Evaluate cost for plot
GlifeND = sum(conf.Acost*iFin(:,1:end-1))*AlpL; % Cost of infected population
GeconND = (1 -(sum(conf.Bcost*(sFin(:,1:end-1) + rFin(:,1:end-1))).*(1-uFin))).*AlpE; % Economic cost
GsocND = (uFin/uM).^2*AlpS; % Cost of lockdown on societal impact

GlifeCum = trapz(tFin(1:end-1),GlifeND)/Tmax; % Cumulative cost of infected population
GeconCum = trapz(tFin(1:end-1),GeconND)/Tmax; % Cumulative economic cost
GsocCum = trapz(tFin(1:end-1),GsocND)/Tmax; % Cumulative cost of lockdown on societal impact

daysLostCum = 100*GeconCum/AlpE; % Number of days of lost economic activity

mortality = [IFRy 0; 0 IFRo];
lifeLost = sum(mortality*(sFin(:,1) - sFin(:,end))); % Total mortality

TsaturHosp = 0*tFin; TsaturHosp(CostriFin>0.95) = 1;
TsaturHospDys = round(max(tFin(TsaturHosp==1)) - min(tFin(TsaturHosp==1))); % Total number of days of Hospital saturation

Uopt = uFin;
Tuopt = tFin(1:end-1); 
%% Plotting 
% close all; clc
np = 3;
mp = 2;
xlima = [0,Tmax];

% Open loop dynamics
[S,I,R,t] = sirOL;

% Closed loop and visualisation 
sFinVf = [solution.states.S1.value; solution.states.S2.value];
iFinVf = [solution.states.I1.value; solution.states.I2.value];
[Svf,Ivf,Sdotvf,Idotvf] = sirVf(sFinVf,iFinVf);

pst = 1;
AxthicksFnt = 13;
figure('units','normalized','outerposition',[0 0 .5 1])

subplot(np,mp,1);
plot(t,I(:,1)/Np,'LineWidth',1.5,'Color',reddish);grid on; hold on
plot(t,I(:,2)/Np,'--','LineWidth',1.5,'Color',reddish);grid on
set(gcf,'renderer','Painters')
pbaspect([2 1 1]);
xlabel('Time'); grid on
xlim(xlima)
set(gca,'FontSize', AxthicksFnt);
set(gcf,'color','w');
title('Uncontrolled'); 
legend('$I_y/N$','$I_o/N$','Interpreter','latex')

subplot(np,mp,2);
plot(S(:,1)/Np,I(:,1)/Np,'LineWidth',2,'Color',reddish);hold on
plot(S(:,2)/Np,I(:,2)/Np,'--','LineWidth',2,'Color',reddish);grid on
xlabel('Susceptible');
ylabel('Infected');
pbaspect([2 1 1]);
hold off;
set(gca,'FontSize', AxthicksFnt);
title('Phase Space Uncontrolled, Germany'); 
legend('Young','Old')

subplot(np,mp,3);
yyaxis left 
hold on
plot(tFin,iFin(1,:),'-','LineWidth',1.5,'Color',reddish);
plot(tFin,iFin(2,:),'-.','LineWidth',1.5,'Color',reddish);
xlabel('Time'); grid on
set(gca,'XTickLabel',[]); grid on
yyaxis right  
plot(tFin,CostriFin,'k','MarkerSize',5,'LineWidth',1.5); hold on 
stairs(tFin(1:end-1),uFin/uM,'k','LineWidth',1.5)
pbaspect([2 1 1]);
xlim(xlima)
ylim([0 1])
set(gca,'FontSize', AxthicksFnt); grid on, box on
title('Controlled'); 
legend('$I_y/N$','$I_o/N$','$I_C$','$u/u_M$','Interpreter','latex')
plt = gca;
plt.YAxis(1).Color = reddish;
plt.YAxis(2).Color = 'k';


subplot(np,mp,4);
xBox = [0, 0, 1, 1, 0];
yBox = [0, Icus, Icus, 0, 0];
plot(sFin(1,:)',iFin(1,:)','LineWidth',2,'Color',reddish);hold on; grid on
plot(sFin(2,:)',iFin(2,:)','--','LineWidth',2,'Color',reddish);
hold on
quiver(Svf(1,1:pst:end),Ivf(1,1:pst:end),Sdotvf(1,1:pst:end),Idotvf(1,1:pst:end),'k');
quiver(Svf(2,1:2*pst:end),Ivf(2,1:2*pst:end),Sdotvf(2,1:2*pst:end),Idotvf(2,1:2*pst:end),'k');
xlabel('Susceptible');
ylabel('Infected');
pbaspect([2 1 1]);
hold off;
set(gca,'FontSize', AxthicksFnt);
title('Phase Space Controlled'); 
legend('Young','Old','Uncontrolled SIR','Location','northwest')

subplot(np,mp,5); 
yyaxis left 
plot(tFin(1:end-1),GlifeND,'-.','LineWidth',1.5,'Color',blueish);hold on
plot(tFin(1:end-1),GeconND,'LineWidth',1.5,'Color',blueish);
plot(tFin(1:end-1),GsocND,'LineWidth',1.5,'Color',blueish); 
yyaxis right 
stairs(tFin(1:end-1),uFin/uM,'k','LineWidth',1)
set(gcf,'renderer','Painters')
pbaspect([2 1 1]);
xlabel('Time'); grid on
xlim(xlima)
set(gca,'FontSize', AxthicksFnt);
set(gcf,'color','w');
title('Cost'); 
legend('$G_{\rm{life}}$','$G_{\rm{econ}}$','$G_{\rm{soc}}$','$u/u_M$','Interpreter','latex')
plt = gca;
plt.YAxis(1).Color = blueish;
plt.YAxis(2).Color = 'k';

subplot(np,mp,6);
plot(sFin(2,:)',iFin(2,:)','--','LineWidth',2,'Color',reddish);
hold on
quiver(Svf(2,1:pst:end),Ivf(2,1:pst:end),Sdotvf(2,1:pst:end),Idotvf(2,1:pst:end),0,'k');
xlabel('Susceptible');
ylabel('Infected');
pbaspect([2 1 1]);
hold off; grid on;
set(gca,'FontSize', AxthicksFnt);
title('Controlled Old'); 
legend('Old','Uncontrolled SIR','Location','south')

% Plot parameters 
Data = [Tmax, Py, Po, Icus, Correc, AlpL, AlpE, AlpS];
VarNames = {'tau', 'tmax', 'py', 'po', 'Icus','Correc','alpL','alpE','alpS'};
fprintf(1, '%s\t\t\t%s\t%s\t\t%s\t\t%s\t\t%s\t%s\t%s\t%s\t\n', VarNames{:});
fprintf(1, '%d\t\t%d\t\t%.3f\t%.4f\t%.4f\t\t%.2f\t%.3f\t%.3f\t%.3f\n', Data');

GlifeCum

GsocCum

GeconCum

lifeLost

daysLostCum

TsaturHospDys
