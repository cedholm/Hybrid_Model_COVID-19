function [COVIDParameters,fvalues,ExitFlags, endpoints] = BC_Phase1_Final_MultiStart(NoStartPoints, Tstart, Tend, place, testnumber) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Karen K. L. Hwang, Christina J. Edholm, Omar Saucedo, Linda J. S. Allen, Nika Shakiba

% BC_Phase1_Final_MultiStart.m

% To run this code you need to give how many Multistart runs you want
% NoStartPoints, the day to start Tstart, the day to end Tend, the place
% choose 'NY', 'ND', or 'BC, and what test this is 1,2,3....
%
% This code calls the COVID_Model_setb3g2mu -- ODE equations
% 
% This code will plot results, calculate R0, AIC, and the minimization
% functional value. Along with possible parameters to minimize.
% 
%
% You can change the LowerBounds and UpperBounds below for the parameters
% we are estimating, beta1, beta2, muI2, and p. This can be chnaged!
% 
%
% June 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Upper and Lower bounds for parameters you are fitting. If you want to change 
% the parameter you need to edit here and in the COVID_Model.m -- only parameters in this code are p and inital consitions. 
% Initial conditions for Asymp and Exposed are set in the next section so
% delete these and set then to z values in the two function COVID_RUN_ODE45
% and COVIDODE45_Plotoutputs.

%Parameter vector we are approximating:
% z = [ beta1, beta2, p, d1] 

LowerBounds=[0.001, 0.001, 0.1, 0.2];
UpperBounds = [2.5, 2.5, 0.9, 0.45];


xstart=.5*(LowerBounds+UpperBounds); %start in the middle

%% Set-up the Population and Time - BC, NY, or ND

%Things to set T, time, day, last date time, Intial Asymp and Exposed

bc=strcmp(place,'BC');
ny=strcmp(place,'NY');
it=strcmp(place,'IT');

if it == 1                                          %% Italy Data
Data=readtable('ItalyApr52020.xlsx','Range','C1:F47');
Data=table2array(Data);
CumCases=Data(:,1)';
CumDeaths=Data(:,2)';
TimeSteps=Data(:,3)';
CumRecov=Data(:,4)';
TotalPop=60000000;
InitialExposedSilent=0;
InitialAsymptomaticSilent=0;
InitialRecoveredSilent=0;
InitialExposedSymptomatic=0;
InitialAsymptomaticSymptomatic=0;
InitialRecoveredSymptomatic=0;

elseif bc == 1                                          %% British Columbia Data
Data=readtable('BC200823_WithRecover.xlsx','Range','C1:G184');
Data=table2array(Data);
CumCases=Data(:,1)';
CumDeaths=Data(:,2)';
TimeSteps=Data(:,3)';
DailyRecovered=Data(:,4)';
CumRecov=Data(:,5)';
TotalPop=5110917;
InitialExposedSilent=0;
InitialAsymptomaticSilent=0;
InitialRecoveredSilent=0;
InitialExposedSymptomatic=0;
InitialAsymptomaticSymptomatic=0;
InitialRecoveredSymptomatic=0;

elseif ny == 1                                        %% New York Data
Data=readtable('NYMay302020.xlsx','Range','D1:F92');
Data=table2array(Data);
CumCases=Data(:,1)';
CumDeaths=Data(:,2)';
TimeSteps=Data(:,3)';
TotalPop=19440469;
InitialExposedSilent=100;
InitialAsymptomaticSilent=100;
InitialRecoveredSilent=0;
InitialExposedSymptomatic=100;
InitialAsymptomaticSymptomatic=100;
InitialRecoveredSymptomatic=0;

else                                                  %% North Dakota Data
Data=readtable('NDMay302020.xlsx','Range','D1:F82');
Data=table2array(Data);
CumCases=Data(:,1)';
CumDeaths=Data(:,2)';
TimeSteps=Data(:,3)';
TotalPop=761723;
InitialExposedSilent=100;
InitialAsymptomaticSilent=100;
InitialRecoveredSilent=0;
InitialExposedSymptomatic=100;
InitialAsymptomaticSymptomatic=100;
InitialRecoveredSymptomatic=0;
end

tspan = Tstart:1:Tend;

datastart=find(TimeSteps==Tstart);

dataend=find(TimeSteps==Tend);

TimeStepsfit=TimeSteps(datastart:dataend);
CumCasesfit=CumCases(datastart:dataend);
CumDeathsfit=CumDeaths(datastart:dataend);
CumRecovfit=CumRecov(datastart:dataend);

InitialRecoveredSymptomatic=CumRecov(datastart);

if Tstart == 1
    InitialInfectedSymptomatic=CumCases(datastart);
else
    InitialInfectedSymptomatic=CumCases(datastart+1)-CumCases(datastart);
end
%% MultiStart and fmincon - Fitting Part - Parallelization

problem = createOptimProblem('fmincon','objective',@COVID_RUN_ODE45...
,'x0',xstart,'lb',LowerBounds,'ub',UpperBounds);%,'Aineq',A,'bineq',b)%,'Aeq',Aeq,'beq',beq);

problem.options = optimoptions(problem.options,'MaxFunEvals',9999,'MaxIter',9999);%,'TolFun',0,'TolCon',0)
%problem.options = optimoptions(problem.options,'MaxFunEvals',inf,'MaxIter',inf,'TolFun',1e-10,'TolCon',0,'TolX',0,'MaxFunEvals',999999)

numstartpoints=NoStartPoints; 

% %  ms=MultiStart('Display','iter');    %defines a multistart problem

ms=MultiStart('UseParallel',true,'Display','iter');       %defines a parallel multistart problem

%parpool %accesses the cores for parallel on your computer (laptop goes for 2-8, can be more specific)

[b,fval,exitflag,output,manymins]=run(ms,problem,numstartpoints);  %runs the multistart 

% the following takes solutions from manymins and makes a matrix out of them


for i=1:length(manymins)
    COVIDParameters(i,:)=manymins(i).X;
end

for i=1:length(manymins)
    fvalues(i)=manymins(i).Fval;
end

for i=1:length(manymins)
    ExitFlags(i)=manymins(i).Exitflag;
end
ExitFlags=ExitFlags';
fvalues=fvalues';

beep on;
beep
beep
beep
beep
beep
beep

%delete(gcp('nocreate'))  %turns off the parallel feature



% Output COVIDParameters, fvalues, and exitflags into an excel file
header={'Exit Flags','fvalues','beta1','beta2','p','d1'};
filename =['COVID_ParamEst_',place,'_TestNumber_',num2str(testnumber),'.xlsx'];
writecell(header,filename,'Sheet',1,'Range','A1');
writematrix(ExitFlags,filename,'Sheet',1,'Range','A2');
writematrix(fvalues,filename,'Sheet',1,'Range','B2');
writematrix(COVIDParameters,filename,'Sheet',1,'Range','C2');



%% Plot the "best" solution - Calculate R0 and AIC

%%Outputs state variables for "best" fit, also gives R0, AIC, and
%%objective function value. 
[yy]=COVIDODE45_Plotoutputs(COVIDParameters(1,:)); 

S1p = yy(:,1); E1p = yy(:,2); A1p = yy(:,3); R1p = yy(:,4); S2p = yy(:,5); E2p = yy(:,6); A2p = yy(:,7); I2p = yy(:,8); R2p = yy(:,9); CIp=yy(:,10); CDp=yy(:,11);


endpoints=[S1p(length(tspan)),E1p(length(tspan)),A1p(length(tspan)), R1p(length(tspan)),S2p(length(tspan)),E2p(length(tspan)),A2p(length(tspan)), I2p(length(tspan)), R2p(length(tspan)),CIp(length(tspan)),CDp(length(tspan))];


% Output endpoints into an excel file
header={'S1','E1','A1', 'R1', 'S2','E2','A2', 'I2', 'R2','CI','CD'};
filename1 =['COVID_ParamEst_endpoints_',place,'_TestNumber_',num2str(testnumber),'.xlsx'];
writecell(header,filename1,'Sheet',1,'Range','A1');
writematrix(endpoints,filename1,'Sheet',1,'Range','A2');

% Create Plots

        figure(1)
        tiledlayout(1,3)
        nexttile
        hold all
        plot(tspan,CIp)
        scatter(TimeStepsfit,CumCasesfit,10,'filled')
        title('Cumulative Cases')
        xlabel('Days')
        nexttile
        hold all
        plot(tspan,CDp)
        scatter(TimeStepsfit,CumDeathsfit,10,'filled')
        title('Cumulative Deaths')
        xlabel('Days')
        savefig(['DataCompare_',place,'_TestNumber_',num2str(testnumber),'.fig'])
        nexttile
        hold all
        plot(tspan,R2p)
        scatter(TimeStepsfit,CumRecovfit,10,'filled')
        title('Cumulative Recovered')
        xlabel('Days')
        savefig(['DataCompare_',place,'_TestNumber_',num2str(testnumber),'.fig'])
        
        figure(2)
        tiledlayout(2,5)
        nexttile
        plot(tspan,S1p)
        title('SilentSpreader S1')
        xlabel('Days')
        nexttile
        plot(tspan,E1p)
        title('SilentSpreader E1')
        xlabel('Days')
        nexttile
        plot(tspan,A1p)
        title('SilentSpreader A1')
        xlabel('Days')
        nexttile
        plot(tspan,R1p)
        title('SilentSpreader R1')
        xlabel('Days')
        nexttile
        nexttile
        plot(tspan,S2p)
        title('SymptomaticSpreader S2')
        xlabel('Days')
        nexttile
        plot(tspan,E2p)
        title('SymptomaticSpreader E2')
        xlabel('Days')
        nexttile
        plot(tspan,A2p)
        title('SymptomaticSpreader A2')
        xlabel('Days')
        nexttile
        plot(tspan,I2p)
        title('SymptomaticSpreader I2')
        xlabel('Days')
        nexttile
        plot(tspan,R2p)
        title('SymptomaticSpreader R2')
        xlabel('Days')
        savefig(['AlltheStateVariables',place,'_TestNumber_',num2str(testnumber),'.fig'])
      
        
function value=COVID_RUN_ODE45(z)

p=z(3);                                                 %fraction of the population that are SilentSpreader

%% Initial Conditions -- Based on T and place


R10=InitialRecoveredSilent;                             %Susceptible individuals for SilentSpreader
E10=InitialExposedSilent;                               %Exposed individuals but not infectious for SilentSpreader
A10=InitialAsymptomaticSilent;                          %Asymptomatic individuals but infectious for SilentSpreader

R20=InitialRecoveredSymptomatic;                        %Susceptible individuals for SymptomaticSpreaders
E20=InitialExposedSymptomatic;                          %Exposed individuals but not infectious for SymptomaticSpreaders
I20=InitialInfectedSymptomatic;                         %Infected individuals for SymptomaticSpreaders
A20=InitialAsymptomaticSymptomatic;                     %Asymptomatic individuals but infectious for SymptomaticSpreaders

S10=p*(TotalPop-R10-E10-A10-R20-E20-A20-I20);                                         %SilentSpreaders
S20=(1-p)*(TotalPop-R10-E10-A10-R20-E20-A20-I20);                                     %SymptomaticSpreaders

%S10=N10-R10-E10-A10;                                    %Removed individuals for SilentSpreader
%S20=N20-R20-E20-A20-I20;                                %Removed individuals for SymptomaticSpreaders

N10=S10+E10+A10+R10;
N20=S20+E20+A20+I20+R20;

CI0=CumCases(datastart);
D0=CumDeaths(datastart);

initialvalues = [S10,E10,A10,R10,S20,E20,A20,I20,R20,CI0,D0]; 

%% ODE45 Solver

[t,y] = ode45(@(t,y) COVID_Model_setb3g2mu_BC(t,y,z),tspan,initialvalues);

      R=y(:,9); CI=y(:,10); CD=y(:,11);
       
%% Data to Fit: TimeSteps recorded days

CIts=zeros(1,length(TimeStepsfit));                 %Cumulative Cases at data points
Dts=zeros(1,length(TimeStepsfit));                  %Cumulative Deaths at data points
RIts= zeros(1,length(TimeStepsfit));                %Recovered per day at data points

for kk=1:length(TimeStepsfit)                       %pulling model values for TimeSteps
    CIts(kk)=CI(TimeStepsfit(kk)-(Tstart-1));
    Dts(kk)=CD(TimeStepsfit(kk)-(Tstart-1));
    RIts(kk)=R(TimeStepsfit(kk)-(Tstart-1));
end

value = 0;                                          % reset just in case

%Difference Data and Model

diff1 = CIts - CumCasesfit;

diff2 = Dts - CumDeathsfit;

diff3 = RIts - CumRecovfit;


%Objective function which is minimized"
value = norm(diff1,2)/norm(CumCasesfit) + norm(diff2,2)/norm(CumDeathsfit)+ norm(diff3,2)/norm(CumRecovfit);

end


function y = COVIDODE45_Plotoutputs(z)

p=z(3);                                             %fraction of the population that are SilentSpreader

%% Initial Conditions -- Based on T and place

R10=InitialRecoveredSilent;                             %Susceptible individuals for SilentSpreader
E10=InitialExposedSilent;                               %Exposed individuals but not infectious for SilentSpreader
A10=InitialAsymptomaticSilent;                          %Asymptomatic individuals but infectious for SilentSpreader

R20=InitialRecoveredSymptomatic;                        %Susceptible individuals for SymptomaticSpreaders
E20=InitialExposedSymptomatic;                          %Exposed individuals but not infectious for SymptomaticSpreaders
I20=InitialInfectedSymptomatic;                         %Infected individuals for SymptomaticSpreaders
A20=InitialAsymptomaticSymptomatic;                     %Asymptomatic individuals but infectious for SymptomaticSpreaders

S10=p*(TotalPop-R10-E10-A10-R20-E20-A20-I20);                                         %SilentSpreaders
S20=(1-p)*(TotalPop-R10-E10-A10-R20-E20-A20-I20);                                     %SymptomaticSpreaders

%S10=N10-R10-E10-A10;                                    %Removed individuals for SilentSpreader
%S20=N20-R20-E20-A20-I20;                                %Removed individuals for SymptomaticSpreaders

N10=S10+E10+A10+R10;
N20=S20+E20+A20+I20+R20;

CI0=CumCases(datastart);
D0=CumDeaths(datastart);

initialvalues = [S10,E10,A10,R10,S20,E20,A20,I20,R20,CI0,D0]; 

%% ODE45 Solver

[t,y] = ode45(@(t,y) COVID_Model_setb3g2mu_BC(t,y,z),tspan,initialvalues);

       R=y(:,9); CI=y(:,10); CD=y(:,11);
       
%% Data to Fit: TimeSteps recorded days

CIts=zeros(1,length(TimeStepsfit));                 %Cumulative Cases at data points
Dts=zeros(1,length(TimeStepsfit));                  %Cumulative Deaths at data points
RIts= zeros(1,length(TimeStepsfit));                %Recovered per day at data points

for kk=1:length(TimeStepsfit)                       %pulling model values for TimeSteps
    CIts(kk)=CI(TimeStepsfit(kk)-(Tstart-1));
    Dts(kk)=CD(TimeStepsfit(kk)-(Tstart-1));
    RIts(kk)=R(TimeStepsfit(kk)-(Tstart-1));
end

value = 0;                                          % reset just in case

%Difference Data and Model

diff1 = CIts - CumCasesfit;

diff2 = Dts - CumDeathsfit;

diff3 = RIts - CumRecovfit;

%Objective function which is minimized"
value = norm(diff1,2)/norm(CumCasesfit) + norm(diff2,2)/norm(CumDeathsfit)+ norm(diff3,2)/norm(CumRecovfit);

fprintf('Case Minimization Value=%s\n',norm(diff1,2)/norm(CumCasesfit));

fprintf('Death Minimization Value=%s\n',norm(diff2,2)/norm(CumDeathsfit));

fprintf('Recovered Minimization Value=%s\n',norm(diff3,2)/norm(CumRecovfit));

fprintf('Total Minimization Value=%s\n',value);

% Calculate R0
N=N10+N20;
beta1=z(1);
beta2=z(2);
beta3=0.3;
delta1= z(4);
delta2= 1/2.3;
gamma2=0.075;
mu2=0.004;
P1t=S10/N;
P2t=S20/N;
Rtt=P1t*(beta1/delta1)+(P2t)*(beta2/delta2+beta3/(gamma2+mu2));

R0=(beta1*N10)/(delta1*N)+(beta2*N20)/(delta2*N)+(beta3*N20)/((gamma2+mu2)*N);

fprintf('Ro=%s\n',R0);
fprintf('Rtt=%s\n',Rtt);

%AIC

%Stores the fitted parameters. 
Fitted_Parameters=length(z); 

error_in_data_C=norm(diff1,2);
error_in_data_D=norm(diff2,2);

%Computes the AIC 
AIC_Cases=length(CumCasesfit)*(log(error_in_data_C/length(CumCasesfit)))+2*Fitted_Parameters; 
AIC_Deaths=length(CumDeathsfit)*(log(error_in_data_D/length(CumDeathsfit)))+2*Fitted_Parameters; 
AIC_Recov=length(CumRecovfit)*(log(error_in_data_D/length(CumRecovfit)))+2*Fitted_Parameters; 

fprintf('AIC Cases=%s\n',AIC_Cases);

fprintf('AIC Deaths=%s\n',AIC_Deaths);

fprintf('AIC Recovered=%s\n',AIC_Recov);


end


end
