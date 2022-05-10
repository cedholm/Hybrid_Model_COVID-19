% Written by: Karen K. L. Hwang, Christina J. Edholm, Omar Saucedo, Linda J. S. Allen, Nika Shakiba

% Mar 24, 2021 
% The purpose of this code is plot the SEAIR model compartments according to
% the fitted parameters of phase 1 in BC. 

close all 
clear all 

p=0.6754;            %fraction of the population that are SilentSpreader
Data = xlsread('210323 BC_Data_WithRecovered.xlsx',1,'C2:G378'); 
CumCases=[Data(:,3),Data(:,1)]; 
CumDeaths=[Data(:,3),Data(:,2)]; 
CumRecover=[Data(:,3),Data(:,5)]; 

%% Initial Conditions 
TotalPop=5110917;

R10=0;              %Susceptible individuals for SilentSpreader
E10=0;              %Exposed individuals but not infectious for SilentSpreader
A10=0;              %Asymptomatic individuals but infectious for SilentSpreader

R20=0;              %Susceptible individuals for SymptomaticSpreaders
E20=0;              %Exposed individuals but not infectious for SymptomaticSpreaders
I20=1;              %Infected individuals for SymptomaticSpreaders
A20=0;              %Asymptomatic individuals but infectious for SymptomaticSpreaders

S10=p*(TotalPop-R10-E10-A10-R20-E20-A20-I20);            %SilentSpreaders
S20=(1-p)*(TotalPop-R10-E10-A10-R20-E20-A20-I20);        %SymptomaticSpreaders

N10=S10+E10+A10+R10;
N20=S20+E20+A20+I20+R20;

CI0=CumCases(1,1);
D0=CumDeaths(1,1);

initialvalues = [S10,E10,A10,R10,S20,E20,A20,I20,R20,CI0,D0]; 

%% ODE45 Solver
tspan=1:length(Data(:,1)); 
[t,y] = ode45(@(t,y) WAMBCovidModel(t,y),tspan,initialvalues);

R=y(:,9); CI=y(:,10); CD=y(:,11);


ModelOutputDay400=y(377,1:11); 

%% Plots
Days=377;
tData=Data(:,3);
tData=tData(tData<=400); 

figure
plot(tspan(1:Days),CI(1:Days),'-b','LineWidth',3)

hold on 
plot(tData,CumCases(1:length(tData),2),'Marker','.','Color',[1 0 0],...
                'MarkerSize',20,'LineStyle','none')
legend({'Cases', 'CCData'},'FontSize',12)
title({'BC Covid19'},...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',16,...
'FontName','Times')
xlabel('Days',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,...
'FontName','Times')
ylabel({'Cumulative Cases'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',14,...
'FontName','Times')


figure
plot(tspan(1:Days),CD(1:Days),'-b','LineWidth',3)

hold on 
plot(tData,CumDeaths(1:length(tData),2),'Marker','.','Color',[1 0 0],...
                'MarkerSize',20,'LineStyle','none')
legend({'Deaths', 'CData'},'FontSize',12)
title({'BC Covid19'},...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',16,...
'FontName','Times')
xlabel('Days',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',14,...
'FontName','Times')
ylabel({'Cumulative Deaths'},...
'FontUnits','points',...
'interpreter','latex',...
'FontSize',14,...
'FontName','Times')

figure
plot(tspan(1:Days),R(1:Days),'-b','LineWidth',3)

hold on 
plot(tData,CumRecover(1:length(tData),2),'Marker','.','Color',[1 0 0],...
                'MarkerSize',20,'LineStyle','none')
legend({'Recovered', 'RData'},'FontSize',12)
title({'BC Covid19'},...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',16,...
'FontName','Times')
xlabel(''-)
ylabel({'Cumulative Recovered'})


