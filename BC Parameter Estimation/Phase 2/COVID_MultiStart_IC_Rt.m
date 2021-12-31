function [COVIDPaRameters,fvalues,ExitFlags, endpoints] = COVID_MultiStart_IC_Rt(NoStartPoints, Tstart, Tend, testnumber) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COVID_MultiStart_IC_Rt.m
%Christina Edholm edited by Karen Hwang
%
% To run this code you need to give how many Multistart runs you want
% NoStartPoints, the day to start Tstart, the day to end Tend
%
% Manually enter the Asymp IC
%
%
% This code calls the COVID_Model_IC_Rt -- ODE equations
% 
% This code will plot results, calculate R0, AIC, Rts and the minimization
% functional value. Along with possible parameters to minimize. (Working on
% AIC and R0)
% 
%
% You can change the LowerBounds and UpperBounds below for the parameters
% we are estimating, beta1, beta2, muI2, and p. This can be changed!
% 
%
% October 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Upper and Lower bounds for parameters you are fitting. If you want to change 
% the parameter you need to edit here and in the COVID_Model.m -- only parameters in this code are p and inital consitions. 
% Initial conditions for Asymp and Exposed are set in the next section so
% delete these and set then to z values in the two function COVID_RUN_ODE45
% and COVIDODE45_Plotoutputs.

%Parameter vector we are approximating:
% z = [ beta1, beta2, beta3, p] 

LowerBounds=[0.3, 0.9, 0.0003, 0.9];
UpperBounds = [0.35, 1, 0.0005, 0.95];											

%xstart=.5*(LowerBounds+UpperBounds); %start in the middle
xstart=[0.305815249,0.983315516,0.000386189,0.937527221];


%% Set-up the Population and Time - BC

%Things to set T, time, day, last date time, Intial Asymp and Exposed

%% British Columbia Data

Data=readtable('BC200823_WithRecover.xlsx','Range','C1:G184');
Data=table2array(Data);
CumCases=Data(:,1)';
CumDeaths=Data(:,2)';
TimeSteps=Data(:,3)';
DailyRecovered=Data(:,4)';
CumRecov=Data(:,5)';
TotalPop=5110917;
P2TotalPop = 5110152.8; %Sum of endpoint SEAIR From phase 1
InitialExposedSilent=1495.101121;
% InitialAsymptomaticSilent=942.8193927;
InitialRecoveredSilent=2408.288725;
InitialExposedSymptomatic=718.3242468;
% InitialAsymptomaticSymptomatic=391.6875341;
InitialInfectedSymptomatic=778.2775254;
InitialRecoveredSymptomatic=418.7481318;

tspan = Tstart:1:Tend;

datastart=find(TimeSteps==Tstart);

dataend=find(TimeSteps==Tend);

TimeStepsfit=TimeSteps(datastart:dataend);
CumCasesfit=CumCases(datastart:dataend);
CumDeathsfit=CumDeaths(datastart:dataend);
CumRecovfit=CumRecov(datastart:dataend);

InitialRecoveredSymptomatic=CumRecov(datastart);

count=1;

%% HERE CHANGE THE INTIAL CONDITIONS

for InitialAsymptomaticSilent=[942.8193927];
    
   for InitialAsymptomaticSymptomatic=[391.6875341];

%% MultiStart and fmincon - Fitting Part - Parallelization - not many comments ask Christina for clarification if you want.

[COVIDParameters,ExitFlags,fvalues] = RunMultiFmincon(xstart,LowerBounds,UpperBounds)

%delete(gcp('nocreate'))  %turns off the parallel feature

[Rout]=CalculateRt(COVIDParameters);

% Output COVIDParameters, fvalues, and exitflags into an excel file
header={'InitialAsymptomaticSilent','InitialAsymptomaticSymptomatic'};
header2={'Exit Flags','fvalues','beta1','beta2','beta3','p','min(Rt)','min(Rtt)','min(Rttt)','min(Rtttt)','median(Rt)','median(Rtt)','median(Rttt)','median(Rtttt)', 'max(Rt)','max(Rtt)','max(Rttt)','max(Rtttt)'};
filename =['COVID_ParamEst_BC_TestNumber_',num2str(testnumber),'.xlsx'];
writecell(header,filename,'Sheet',count,'Range','B1');
writematrix([InitialAsymptomaticSilent,InitialAsymptomaticSymptomatic],filename,'Sheet',count,'Range','B2');

writecell(header2,filename,'Sheet',count,'Range','A3');
writematrix(ExitFlags,filename,'Sheet',count,'Range','A4');
writematrix(fvalues,filename,'Sheet',count,'Range','B4');
writematrix(COVIDParameters,filename,'Sheet',count,'Range','C4');
writematrix(Rout,filename,'Sheet',count,'Range','H4');



%% Plot the "best" solution - Calculate R0 and AIC

%%Outputs state variables for "best" fit, also gives R0, AIC, and
%%objective function value. 
[yy]=COVIDODE45_Plotoutputs(COVIDParameters(1,:)); 

S1p = yy(:,1); E1p = yy(:,2); A1p = yy(:,3); R1p = yy(:,4); S2p = yy(:,5); E2p = yy(:,6); A2p = yy(:,7); I2p = yy(:,8); R2p = yy(:,9); CIp=yy(:,10); CDp=yy(:,11);

endpoints=[S1p(length(tspan)),E1p(length(tspan)),A1p(length(tspan)), R1p(length(tspan)),S2p(length(tspan)),E2p(length(tspan)),A2p(length(tspan)), I2p(length(tspan)), R2p(length(tspan)),CIp(length(tspan)),CDp(length(tspan))];


% Output endpoints into an excel file
header3={'InitialAsymptomaticSilent','InitialAsymptomaticSymptomatic'};
header4={'S1','E1','A1', 'R1', 'S2','E2','A2', 'I2', 'R2','CI','CD'};
filename1 =['COVID_ParamEst_endpoints_BC_TestNumber_',num2str(testnumber),'.xlsx'];
writecell(header3,filename1,'Sheet',count,'Range','B1');
writematrix([InitialAsymptomaticSilent,InitialAsymptomaticSymptomatic],filename1,'Sheet',count,'Range','B2');
writecell(header4,filename1,'Sheet',count,'Range','A3');
writematrix(endpoints,filename1,'Sheet',count,'Range','A4');

% Create Plots

        figure(count)
        tiledlayout(3,5)
        nexttile
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
        nexttile
        hold all
        plot(tspan,R2p)
        scatter(TimeStepsfit,CumRecovfit,10,'filled')
        title('Cumulative Recovered')
        xlabel('Days')
        nexttile
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
        savefig(['Everything_BC_SweepNumber_',num2str(count),'_TestNumber_',num2str(testnumber),'.fig'])
        
 count=count+1;       
    end
end

        
%% Functions       
   

    function [COVIDParameters,ExitFlags,fvalues] = RunMultiFmincon(xstart,LowerBounds,UpperBounds)
        
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

end




function value=COVID_RUN_ODE45(z)

p=z(4);                                                 %fraction of the population that are SilentSpreader

%% Initial Conditions -- Based on T and place


R10=InitialRecoveredSilent;                             %Susceptible individuals for SilentSpreader
E10=InitialExposedSilent;                               %Exposed individuals but not infectious for SilentSpreader
A10=InitialAsymptomaticSilent;                          %Asymptomatic individuals but infectious for SilentSpreader

R20=InitialRecoveredSymptomatic;                        %Susceptible individuals for SymptomaticSpreaders
E20=InitialExposedSymptomatic;                          %Exposed individuals but not infectious for SymptomaticSpreaders
I20=InitialInfectedSymptomatic;                         %Infected individuals for SymptomaticSpreaders
A20=InitialAsymptomaticSymptomatic;                     %Asymptomatic individuals but infectious for SymptomaticSpreaders

S10=p*(P2TotalPop-R10-E10-A10-R20-E20-A20-I20);                                         %SilentSpreaders
S20=(1-p)*(P2TotalPop-R10-E10-A10-R20-E20-A20-I20);                                     %SymptomaticSpreaders
%S10=3447000;
%S20=1656000;  
S10+S20

%S10=N10-R10-E10-A10;                                    %Removed individuals for SilentSpreader
%S20=N20-R20-E20-A20-I20;                                %Removed individuals for SymptomaticSpreaders

N10=S10+E10+A10+R10;
N20=S20+E20+A20+I20+R20;

CI0=CumCases(datastart);
D0=CumDeaths(datastart);

initialvalues = [S10,E10,A10,R10,S20,E20,A20,I20,R20,CI0,D0]; 

%% ODE45 Solver

[t,y] = ode45(@(t,y) COVID_Model_IC_Rt(t,y,z),tspan,initialvalues);

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
%value = norm(diff1,2)/norm(CumCasesfit) + norm(diff2,2)/norm(CumDeathsfit)+ norm(diff3,2)/norm(CumRecovfit);

value = norm(diff2,2)/norm(CumDeathsfit);

vec1 = diff1(5:length(diff1));

if sum(vec1<0) > 0
    value = 1000000000;
end

end


function y = COVIDODE45_Plotoutputs(z)

p=z(4);                                             %fraction of the population that are SilentSpreader

%% Initial Conditions -- Based on T and place

R10=InitialRecoveredSilent;                             %Susceptible individuals for SilentSpreader
E10=InitialExposedSilent;                               %Exposed individuals but not infectious for SilentSpreader
A10=InitialAsymptomaticSilent;                          %Asymptomatic individuals but infectious for SilentSpreader

R20=InitialRecoveredSymptomatic;                        %Susceptible individuals for SymptomaticSpreaders
E20=InitialExposedSymptomatic;                          %Exposed individuals but not infectious for SymptomaticSpreaders
I20=InitialInfectedSymptomatic;                         %Infected individuals for SymptomaticSpreaders
A20=InitialAsymptomaticSymptomatic;                     %Asymptomatic individuals but infectious for SymptomaticSpreaders

S10=p*(P2TotalPop-R10-E10-A10-R20-E20-A20-I20);                                         %SilentSpreaders
S20=(1-p)*(P2TotalPop-R10-E10-A10-R20-E20-A20-I20);                                     %SymptomaticSpreaders
%S10=3447000;
%S20=1656000; 

%S10=N10-R10-E10-A10;                                    %Removed individuals for SilentSpreader
%S20=N20-R20-E20-A20-I20;                                %Removed individuals for SymptomaticSpreaders

N10=S10+E10+A10+R10;
N20=S20+E20+A20+I20+R20;

CI0=CumCases(datastart);
D0=CumDeaths(datastart);

initialvalues = [S10,E10,A10,R10,S20,E20,A20,I20,R20,CI0,D0]; 

%% ODE45 Solver

[t,y] = ode45(@(t,y) COVID_Model_IC_Rt(t,y,z),tspan,initialvalues);

    S1 = y(:,1); E1 = y(:,2); A1 = y(:,3); R1 = y(:,4); S2 = y(:,5); E2 = y(:,6); A2 = y(:,7); I2 = y(:,8); R2 = y(:,9); CI=y(:,10); CD=y(:,11); 
    
       
%% Data to Fit: TimeSteps recorded days

CIts=zeros(1,length(TimeStepsfit));                 %Cumulative Cases at data points
Dts=zeros(1,length(TimeStepsfit));                  %Cumulative Deaths at data points
RIts= zeros(1,length(TimeStepsfit));                %Recovered per day at data points

for kk=1:length(TimeStepsfit)                       %pulling model values for TimeSteps
    CIts(kk)=CI(TimeStepsfit(kk)-(Tstart-1));
    Dts(kk)=CD(TimeStepsfit(kk)-(Tstart-1));
    RIts(kk)=R2(TimeStepsfit(kk)-(Tstart-1));
end

value = 0;                                          % reset just in case

%Difference Data and Model

diff1 = CIts - CumCasesfit;

diff2 = Dts - CumDeathsfit;

diff3 = RIts - CumRecovfit;

%Objective function which is minimized"
%value = norm(diff1,2)/norm(CumCasesfit) + norm(diff2,2)/norm(CumDeathsfit)+ norm(diff3,2)/norm(CumRecovfit);

%just fitting deaths
value = norm(diff2,2)/norm(CumDeathsfit);

vec1 = diff1(5:length(diff1));

if sum(vec1<0) > 0
    value = 1000000000;
end

fprintf('Case Minimization Value=%s\n',norm(diff1,2)/norm(CumCasesfit));

fprintf('Death Minimization Value=%s\n',norm(diff2,2)/norm(CumDeathsfit));

fprintf('Recovered Minimization Value=%s\n',norm(diff3,2)/norm(CumRecovfit));

fprintf('Total Minimization Value=%s\n',value);

% Calculate R0
N=N10+N20;
beta1=z(1);
beta2=z(2);
beta3=z(3);
delta1= 0.356823803; 
delta2= 1/2.3;
mu2=0.00314;
gamma2=0.075;

R0=(beta1*N10)/(delta1*N)+(beta2*N20)/(delta2*N)+(beta3*N20)/((gamma2+mu2)*N);
Reff= ((S10+S20)/N)*2.755;

fprintf('Ro=%s\n',R0);
fprintf('Reff=%s\n',Reff);


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


function [Rtoutt] = CalculateRt(COVIDParam) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CalculateRt.m for BC
%Christina Edholm

%loop through all param combinations

for ll=1:length(COVIDParam)

    z=COVIDParam(ll,:);
    
%% Parameters

    p=z(4);                                             %fraction of the population that are SilentSpreader
    beta1=z(1);
    beta2=z(2);
    beta3=z(3);
    delta1= 0.356823803; 
    delta2= 1/2.3;
    mu2=0.00314;
    gamma2=0.075; 

%% Initial Conditions -- For BC

R10=InitialRecoveredSilent;                             %Susceptible individuals for SilentSpreader
E10=InitialExposedSilent;                               %Exposed individuals but not infectious for SilentSpreader
A10=InitialAsymptomaticSilent;                          %Asymptomatic individuals but infectious for SilentSpreader

R20=InitialRecoveredSymptomatic;                        %Susceptible individuals for SymptomaticSpreaders
E20=InitialExposedSymptomatic;                          %Exposed individuals but not infectious for SymptomaticSpreaders
I20=InitialInfectedSymptomatic;                         %Infected individuals for SymptomaticSpreaders
A20=InitialAsymptomaticSymptomatic;                     %Asymptomatic individuals but infectious for SymptomaticSpreaders

S10=p*(P2TotalPop-R10-E10-A10-R20-E20-A20-I20);                                         %SilentSpreaders
S20=(1-p)*(P2TotalPop-R10-E10-A10-R20-E20-A20-I20);                                     %SymptomaticSpreaders
%S10=3447000;
%S20=1656000; 

initialvalues = [S10,E10,A10,R10,S20,E20,A20,I20,R20]; 

%% Run ODE
[t,y] = ode45(@(t,y) f(t,y,z),tspan,initialvalues);

S1 = y(:,1); E1 = y(:,2); A1 = y(:,3); R1 = y(:,4); S2 = y(:,5); E2 = y(:,6); A2 = y(:,7); I2 = y(:,8); R2 = y(:,9);

Pt=(S1./(S1+E1+A1+R1+S2+E2+A2+I2+R2));
Rt=Pt.*(beta1/delta1)+(1-Pt).*(beta2/delta2+beta3/(gamma2+mu2));

Ptt=((S1+E1+A1+R1)./(S1+E1+A1+R1+S2+E2+A2+I2+R2));
Rtt=Ptt.*(beta1/delta1)+(1-Ptt).*(beta2/delta2+beta3/(gamma2+mu2));   

P1t=S1./(S1+E1+A1+R1+S2+E2+A2+I2+R2);
P2t=S2./(S1+E1+A1+R1+S2+E2+A2+I2+R2);
Rttt=P1t.*(beta1/delta1)+(P2t).*(beta2/delta2+beta3/(gamma2+mu2));

P1tt=S1./(P2TotalPop);
P2tt=S2./(P2TotalPop);
Rtttt=P1tt.*(beta1/delta1)+(P2tt).*(beta2/delta2+beta3/(gamma2+mu2));

Rtouttmedian(ll,:)=[median(Rt), median(Rtt), median(Rttt), median(Rtttt)];
Rtouttmin(ll,:)=[min(Rt), min(Rtt), min(Rttt), min(Rtttt)];
Rtouttmax(ll,:)=[max(Rt), max(Rtt), max(Rttt), max(Rtttt)];

Rtoutt=[Rtouttmin,Rtouttmedian,Rtouttmax];

end

end

function dydt = f(t,y,z)

%Parameters

%Parameter Values - COVID

b1=z(1);            %transmission rate for SilentSpreader Asymptomatics
b2=z(2);            %transmission rate for SymptomaticSpreader Asymptomatics
b3=z(3);              %transmission rate for SymptomaticSpreader Infectious
aa1=1/(5.5-2.3);    %rate of transition from exposed to asymptomatic stage for SilentSpreader - COVID -- double check
aa2=1/(5.5-2.3);    %rate of transition from exposed to asymptomatic stage for SymptomaticSpreader - COVID -- double check
d1= 0.356823803;          %rate of transition from asymptomatic to infected stage for SilentSpreader - COVID
d2= 1/2.3;          %rate of transition from asymptomatic to infected stage for SymptomaticSpreader - COVID
mI2=0.00314;           %disease-induced mortality rate for infected SymptomaticSpreader - MERS -- region dependent
g2=0.075;            %removal rate for SymptomaticSpreader - COVID


%leave these alone for now
ep1=1;              %scaling factor for SilentSpreader Asymptomatics
ep2=1;              %scaling factor for SymptomaticSpreader Asymptomatics
rh2=1;              %scaling factor for SymptomaticSpreader Infectious


s1=y(1); e1=y(2); a1=y(3); r1=y(4);             %SilentSpreader classes
s2=y(5); e2=y(6); a2=y(7); y2=y(8); r2=y(9);    %Sypmtomatic classes

n1=s1+e1+a1+r1; n2=s2+e2+a2+y2+r2; %Population sized for all groups

prop=s1/(n1+n2);
prop2=s2/(n1+n2);

s1dot= -1*prop*b1*(ep1*a1) - prop*b2*(ep2*a2)- prop*b3*(rh2*y2);
s2dot= -1*prop2*b1*(ep1*a1) - prop2*b2*(ep2*a2)- prop2*b3*(rh2*y2);
e1dot= prop*b1*(ep1*a1) + prop*b2*(ep2*a2) + prop*b3*(rh2*y2) - aa1*e1;
e2dot= prop2*b1*(ep1*a1) + prop2*b2*(ep2*a2) + prop2*b3*(rh2*y2) - aa2*e2;
a1dot= aa1*e1 - d1*a1;
a2dot= aa2*e2 - d2*a2;
y2dot= d2*a2 - g2*y2 - mI2*y2;
r1dot= d1*a1;
r2dot= g2*y2;

dydt=[s1dot; e1dot; a1dot; r1dot; s2dot; e2dot; a2dot; y2dot; r2dot];

%dydt=dydt'

end





end