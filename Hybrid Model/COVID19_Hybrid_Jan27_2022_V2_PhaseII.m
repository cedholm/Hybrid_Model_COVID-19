% Written by: Karen K. L. Hwang, Christina J. Edholm, Omar Saucedo, Linda J. S. Allen, Nika Shakiba

%Hybrid model combining:
%Time Inhomogeneous Markov chain (Non-homogenous Process = NHP)
%Stochastic Differential Equation (SDNE)

%Random Variables are Discrete
%Monte Carlo Simulation Demographic
%SEIR COVID19 Model
%Computes probability of Extinction

clear all
clc;
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENVIRONMENTAL VARIABILITY PARAMETER VALUES

r_value = 0;
CV_value = 0;
InitialCondition = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARMATERS ESTIMATED USING ODE MODEL-FITTING OR CLINICAL MEASUREMENTS - COVID-19

% Fit parameters
global fraction_silent fraction_symp b1_0 b2_0 b3_0 mI2 aa1 aa2 dd1 dd2 g2 POPULATION N_0 N1_0 N2_0 Rt_0 r1 r2 r3 sigma1 sigma2 sigma3 time dt


% Phase 2 Parameters
fraction_silent = 0.937944626;          %fraction of individuals that are silent spreaders
fraction_symp = 1-fraction_silent;      %fraction of individuals that are symptomatic spreaders
b1_0 = 0.305688791;                     %transmission rate for SilentSpreader Asymptomatic
b2_0 = 0.994680741;                     %transmission rate for SymptomaticSpreader Asymptomatic (Presymptomatic)
b3_0 = 0.000361678;                     %transmission rate for SymptomaticSpreader Infectious
mI2 = 0.00314;                          %disease-induced mortality rate for infected SymptomaticSpreader - COVID-19 -- region dependent


% Measured parameters
aa1 = 1/(5.5-2.3);                      %rate of transition from exposed to asymptomatic stage for SilentSpreader - COVID
aa2 = 1/(5.5-2.3);                      %rate of transition from exposed to asymptomatic stage for SymptomaticSpreader - COVID
dd1 = 0.356907003;                      %rate of transition from asymptomatic to infected stage for SilentSpreader - COVID
dd2 = 1/2.3;                            %rate of transition from asymptomatic to infected stage for SymptomaticSpreader - COVID
g2 = 0.075;                             %removal rate for SymptomaticSpreader - COVID

% POPULATION = 5110917;                 %population size of jurisdiction of interest

% Phase 2 Population
POPULATION = 5110894.661;


% Set up the population
N_0 = round(POPULATION/1);                      % Simulate for a fraction of total population
N1_0 = round(fraction_silent*N_0);              % Silent spreaders at t0
N2_0 = round((1-fraction_silent)*N_0);          % Symptomatic spreaders at t0

% For ICU bed availability, need fraction symptomatic (I2) who end up in
% ICU
frac_ICU = 0.0132;
BC_ICU_beds = 313;

%Calculate Rt at the start of Phase 2
Rt_0 = b1_0*(N1_0)/(dd1*(N_0))+b2_0*(N2_0)/(dd2*(N_0))+b3_0*(N2_0)/((g2+mI2)*(N_0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS
if InitialCondition == 1
    E1_0 = 1496.700603;
    A1_0 = 941.807072;
    R1_0 = 2409.030412;
        
    E2_0 = 719.3204259;
    A2_0 = 391.1772991;
    I2_0 = 779.0445821;
    R2_0 = 418.864848;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME PARAMETERS AND NUMBER OF SIMULATIONS
nsim = 2500;
 
time = 400;           
dt = .0005;            %%%Prob<1, if Check NOT equal zero, then make smaller time steps.
sdt = sqrt(dt);
ts = time/dt;

outbreak = 50;

% Threshold number of E+A+I individuals that requires use of SDE
threshold_hybrid_NHPtoSDE = 100;
threshold_hybrid_SDEtoNHP = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE OUTPUT FILES
% Delete files for output data if they exist

clearvars r_value_print
if mod(r_value,1)
    temp = num2str(r_value);
    for i = 1:length(temp)
        if temp(i) == '.'
            r_value_print(i) = 'p';
        else
            r_value_print(i) = temp(i);
        end
    end
else
    r_value_print = num2str(r_value);
end

clearvars CV_value_print
if mod(CV_value,1)
    temp = num2str(CV_value);
    for i = 1:length(temp)
        if temp(i) == '.'
            CV_value_print(i) = 'p';
        else
            CV_value_print(i) = temp(i);
        end
    end
else
    CV_value_print = num2str(CV_value);
end

outbreak_value_print = num2str(outbreak);

filename1 = ['Hybrid_COVID19_PeakInf_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
if exist(filename1, 'file')==2
  delete(filename1);
end
filename2 = ['Hybrid_COVID19_TimePeak_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
if exist(filename2, 'file')==2
  delete(filename2);
end
filename3 = ['Hybrid_COVID19_NumDeaths_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
if exist(filename3, 'file')==2
  delete(filename3);
end
filename4 = ['Hybrid_COVID19_ProbOutbreak_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
if exist(filename4, 'file')==2
  delete(filename4);
end
% filename5 = ['Hybrid_COVID19_ProbExtinct_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
% if exist(filename5, 'file')==2
%   delete(filename5);
% end
filename6 = ['Hybrid_COVID19_NumInf_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
if exist(filename6, 'file')==2
  delete(filename6);
end
filename7 = ['Hybrid_COVID19_TimeOutbreak_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
if exist(filename7, 'file')==2
  delete(filename7);
end
filename8 = ['Hybrid_COVID19_NumSilentInf_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
if exist(filename8, 'file')==2
  delete(filename8);
end
filename9 = ['Hybrid_COVID19_ProbICUoverload_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
if exist(filename9, 'file')==2
  delete(filename9);
end
filename10 = ['Hybrid_COVID19_NHPerrorfreq_OutbreakThresh_',outbreak_value_print,'_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition),'_PhaseII','.txt'];
if exist(filename10, 'file')==2
  delete(filename10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hybrid ALGORITHM

% Noisiness associated with beta1, 2 and 3
r1 = r_value;
r2 = r1;
r3 = r1;
sigval1 = 2*b1_0*CV_value^2;  
sigval2 = 2*b2_0*CV_value^2; 
sigval3 = 2*b3_0*CV_value^2; 
sigma1 = sqrt(sigval1*r1);              % std in beta1 SDE
sigma2 = sqrt(sigval2*r2);              % std in beta2 SDE
sigma3 = sqrt(sigval3*r3);              % std in beta3 SDE

Error_check=zeros(nsim,time/dt+1);
Error_run=0;

clearvars Time_Outbreak Time_Ext Time_Peak_Infections Num_Deaths Num_Inf Num_SilInf

numOutbreak = 0;
numExtinct = 0;
numICUoverload = 0;

for kk=1:nsim
    kk
    clear t beta1 beta2 beta3 S1 S2 E1 E2 A1 A2 I2 R1 R2

    outbreak_occured = false;
    ICU_overload = false;
    timestep = 1;
    time_curr = 0;
    numDeaths = 0;
    Switch = 0;

    beta1(timestep) = b1_0;
    beta2(timestep) = b2_0;
    beta3(timestep) = b3_0;
    
    E1(timestep)=E1_0;                  E2(timestep)=E2_0;
    A1(timestep)=A1_0;                  A2(timestep)=A2_0;
                                        I2(timestep)=I2_0;
    R1(timestep)=R1_0;                  R2(timestep)=R2_0;
    S1(timestep)=N1_0-E1_0-A1_0-R1_0;   S2(timestep)=N2_0-E2_0-A2_0-I2_0-R2_0;
    
    
    while time_curr <= time && E1(timestep)+E2(timestep)+A1(timestep)+A2(timestep)+I2(timestep) > 0
        if Switch == 0
            % Call NHP

            [timestep,time_curr,beta1,beta2,beta3,S1,S2,E1,E2,A1,A2,I2,R1,R2,numDeaths,Error_check,Switch] =...
            COVID19_Hybrid_NHP_Feb18_2021_V2(timestep,time_curr,kk,beta1,beta2,beta3,S1,S2,E1,E2,A1,A2,I2,R1,R2,numDeaths,Error_check,threshold_hybrid_NHPtoSDE,Switch);
        
        elseif Switch == 1
            %call SDE

            [timestep,time_curr,beta1,beta2,beta3,S1,S2,E1,E2,A1,A2,I2,R1,R2,numDeaths,Switch] =...
            COVID19_Hybrid_SDE_Feb18_2021_V2(timestep,time_curr,beta1,beta2,beta3,S1,S2,E1,E2,A1,A2,I2,R1,R2,numDeaths,threshold_hybrid_SDEtoNHP,Switch);
        end
    end
    
    % Count number of times ICU beds ran out
    for ii = 1:timestep
        if I2(ii)*frac_ICU > BC_ICU_beds && ICU_overload == false
            numICUoverload = numICUoverload + 1;
            ICU_overload = true;
        end
    end
    
    % Record time of outbreak if outbreak occured
    for ii = 1:timestep
        if I2(ii) >= outbreak && outbreak_occured == false
            Time_Outbreak(numOutbreak+1) = dt*ii;
            outbreak_occured = true;
        end

        if E1(ii)+E2(ii)+A1(ii)+A2(ii)+I2(ii) == 0
            [minInf time_minInf] = min(E1+E2+A1+A2+I2);
            Time_Ext(numExtinct+1) = time_minInf*dt;

            numExtinct = numExtinct + 1;
        end
    end

    if outbreak_occured == true
        [peak peaktime] = max(I2);
        Peak_Infections(numOutbreak+1) = peak;
        Time_Peak_Infections(numOutbreak+1) = peaktime*dt;
        Num_Deaths(numOutbreak+1) = numDeaths;
        Num_Inf(numOutbreak+1) = N_0 - S1(timestep) - S2(timestep);
        Num_SilInf(numOutbreak+1) = A1(timestep)+E1(timestep)+R1(timestep);

        numOutbreak = numOutbreak + 1;
    end
    
    if sum(Error_check(kk,:)) > 0
        Error_run = Error_run + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE DATA
dlmwrite(filename1,Peak_Infections,'delimiter','\t','precision',10,'-append')
dlmwrite(filename2,Time_Peak_Infections,'delimiter','\t','precision',10,'-append')
dlmwrite(filename3,Num_Deaths,'delimiter','\t','precision',10,'-append')
dlmwrite(filename4,numOutbreak/nsim,'delimiter','\t','precision',10,'-append')
% dlmwrite(filename5,numExtinct/nsim,'delimiter','\t','precision',10,'-append')
dlmwrite(filename6,Num_Inf,'delimiter','\t','precision',10,'-append')
dlmwrite(filename7,Time_Outbreak,'delimiter','\t','precision',10,'-append')
dlmwrite(filename8,Num_SilInf,'delimiter','\t','precision',10,'-append')
dlmwrite(filename9,numICUoverload/nsim,'delimiter','\t','precision',10,'-append')
dlmwrite(filename10,Error_run/nsim,'delimiter','\t','precision',10,'-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINT DATA SUMMARY
disp((sprintf('----------------------------------------------------------------------')))
disp((sprintf('nsim = %d, time = %d, timestep = %f, thresh = %d\nr = %.2f, CV = %.2f\nE1_0 = %d, E2_0 = %d, A1_0 = %d, A2_0 = %d, I2_0 = %d',nsim,time,dt,outbreak,r_value,CV_value,E1_0,E2_0,A1_0,A2_0,I2_0)))
disp((sprintf('----------------------------------------------------------------------')))
disp((sprintf('Check %d',sum(sum(Error_check)))))
disp((sprintf('----------------------------------------------------------------------')))
disp((sprintf('Rt_0 = %.2f',Rt_0)))
disp((sprintf('Probability of outbreak = %.2f',numOutbreak/nsim)))
disp((sprintf('Probability of ICU overload = %.2f',numICUoverload/nsim)))
disp((sprintf('----------------------------------------------------------------------')))
disp((sprintf('PeakInf \t %.2f \t %.2f',mean(Peak_Infections),std(Peak_Infections))))
disp((sprintf('TimePeakInf \t %.2f \t %.2f',mean(Time_Peak_Infections),std(Time_Peak_Infections))))
disp((sprintf('NumDeaths \t %.2f \t %.2f',mean(Num_Deaths),std(Num_Deaths))))
disp((sprintf('NumInf \t\t %.2f \t %.2f',mean(Num_Inf),std(Num_Inf))))
disp((sprintf('TimeOutbreak \t %.2f \t %.2f',mean(Time_Outbreak),std(Time_Outbreak))))
disp((sprintf('----------------------------------------------------------------------')))

