%Hybrid model combining:
%Time Inhomogeneous Markov chain (Non-homogenous Process = NHP)
%Stochastic Differential Equation (SDNE)

%Random Variables are Discrete
%Monte Carlo Simulation Demographic
%SEIR COVID19 Model
%Computes probability of Extinction

%Edited Feb 15 by Karen
%Changed R0 equation to Rttt equation (Variable name still R0)
%Dragged the Rttt equation down after parameter get updated timesteps
%right before the code for NHP <-> SDE 

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
global fraction_silent fraction_symp R0_fit b1_0 b2_0 b3_0 mI2 aa1 aa2 dd1 dd2 g2 POPULATION N_0 N1_0 N2_0 R0 r1 r2 r3 sigma1 sigma2 sigma3 time dt

% Phase 1 Parameters
% fraction_silent = 0.675469407;          %fraction of individuals that are silent spreaders
% fraction_symp = 1-fraction_silent;      %fraction of individuals that are symptomatic spreaders
% R0_fit = 2.755096228;                   %basic reproduction number
% b1_0 = 0.600106951;                     %transmission rate for SilentSpreader Asymptomatic
% b2_0 = 0.519333416;                     %transmission rate for SymptomaticSpreader Asymptomatic (Presymptomatic)
% b3_0 = 0.3;                             %transmission rate for SymptomaticSpreader Infectious
% mI2 = 0.004;                            %disease-induced mortality rate for infected SymptomaticSpreader - COVID-19 -- region dependent

% Phase 2 Parameters
fraction_silent = 0.937944626;          %fraction of individuals that are silent spreaders
fraction_symp = 1-fraction_silent;      %fraction of individuals that are symptomatic spreaders
R0_fit = 2.755096228;                   %basic reproduction number
b1_0 = 0.305688791;                     %transmission rate for SilentSpreader Asymptomatic
b2_0 = 0.994680741;                     %transmission rate for SymptomaticSpreader Asymptomatic (Presymptomatic)
b3_0 = 0.000361678;                             %transmission rate for SymptomaticSpreader Infectious
mI2 = 0.00314;                            %disease-induced mortality rate for infected SymptomaticSpreader - COVID-19 -- region dependent


% Measured parameters
aa1 = 1/(5.5-2.3);                      %rate of transition from exposed to asymptomatic stage for SilentSpreader - COVID
aa2 = 1/(5.5-2.3);                      %rate of transition from exposed to asymptomatic stage for SymptomaticSpreader - COVID
dd1 = 0.356907003;                      %rate of transition from asymptomatic to infected stage for SilentSpreader - COVID
dd2 = 1/2.3;                            %rate of transition from asymptomatic to infected stage for SymptomaticSpreader - COVID
g2 = 0.075;                             %removal rate for SymptomaticSpreader - COVID

% POPULATION = 5110917;                   %population size of jurisdiction of interest

% Phase 2 Population
POPULATION = 5110894.661;


% Set up the population
N_0 = round(POPULATION/1);                      % Simulate for a fraction of total population
N1_0 = round(fraction_silent*N_0);                 % Silent spreaders at t0
N2_0 = round((1-fraction_silent)*N_0);             % Symptomatic spreaders at t0

% For ICU bed availability, need fraction symptomatic (I2) who end up in
% ICU
frac_ICU = 0.0132;
BC_ICU_beds = 313;

%Nika’s original R0 eqn
%R0 = b1_0*(N1_0)/(dd1*(N_0))+b1_0*(N2_0)/(dd2*(N_0))+b2_0*(N2_0)/((g2+mI2)*(N_0));

%From Phase 1 Final Multistart
%R0=(beta1*N10)/(delta1*N)+(beta2*N20)/(delta2*N)+(beta3*N20)/((gamma2+mu2)*N);

% Updated R0 for COVID Phase 1 (Please check!)
R0 = b1_0*(N1_0)/(dd1*(N_0))+b2_0*(N2_0)/(dd2*(N_0))+b3_0*(N2_0)/((g2+mI2)*(N_0));


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

% What is considered an Outbreak? Maybe based on ICU capacity, maybe based
% on what we actually considered outbreak in Wuhan, etc., can look at the
% impact of threshold on how quickly we go from yellow to red lockdown
% zones
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
    
%     % Rttt for phase 2
%     P1t=S1./(S1+E1+A1+R1+S2+E2+A2+I2+R2);
%     P2t=S2./(S1+E1+A1+R1+S2+E2+A2+I2+R2);
%     R0 = P1t.*(b1_0/dd1)+(P2t).*(b2_0/dd2+b3_0/(g2+mI2));
    
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

%     if kk <= 10
%         cmap=colormap(hsv(10));
%         figure(1)
%         subplot(1,2,1)
%         hold on
%         plot(0:dt:(length(I2)-1)*dt,I2,'Linewidth',1.5,'Color',cmap(kk,:));
%         set(gca,'Linewidth',2)
%         axis tight
%         ylabel('Number of infections')
%         xlabel('Time (days)')
%         grid off
%         box on
% 
%         if outbreak_occured == true
%             subplot(1,2,2)
%             hold on
%             plot(0:dt:(length(I2)-1)*dt,beta2,'Linewidth',1.5,'Color',cmap(kk,:));
%             set(gca,'Linewidth',2)
%             axis tight
%             ylabel('\beta_2(t)')
%             xlabel('Time (days)')
%             grid off
%             box on
%         end
%     end
    
    % Plot the error check
%     figure(3)
%     subplot(1,2,1)
%     c = rand(1,3);
%     plot(0:dt:(length(I2)-2)*dt,Error_check(kk,1:length(I2)-1),'Linewidth',1.5,'color',c)
%     set(gca,'Linewidth',2)
%     ylabel('Error check value')
%     xlabel('Time (days')
%     box on
%     hold on
%     
%     if sum(Error_check(kk,:)) > 0
%         subplot(1,2,2)
%         plot(0:dt:(length(I2)-1)*dt,I2(1:timestep-1),'Linewidth',1.5,'Color',c)
%         set(gca,'Linewidth',2)
%         ylabel('Number of infections')
%         xlabel('Time (days')
%         box on
%     end
%     hold on
    
    if sum(Error_check(kk,:)) > 0
        Error_run = Error_run + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE PLOTS
% figure(1)
% subplot(1,2,2)
% hold on
% plot(0:dt:(length(I2)-1)*dt,b2_0*ones(1,length(I2)),'k:','Linewidth',2);
% 
% str = sprintf(['NHP_COVID19_SamplePathNumInf_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition)]);
% set(gcf,'Name',str,'NumberTitle','off')
% saveas(gcf,strcat(str,'.fig'))
% saveas(gcf,strcat(str,'.png'))
% 
% figure(2)
% subplot(6,1,1)
% [h temp] = hist(Peak_Infections,50);
% b = bar(temp,h/nsim);
% set(b,'FaceColor',[51 51 255]/255,'EdgeColor','k','Linewidth',1.5);
% set(gca,'Linewidth',2)
% set(gca,'Linewidth',2)
% axis tight
% str = sprintf('R_0 = %.2f, I2_0=%i (Mean = %.2f, Std = %f)',R0, I2_0, mean(Peak_Infections),std(Peak_Infections));
% title(str);
% ylabel('Frequency')
% xlabel('Peak infections')
% grid on
% grid minor
% grid minor
% 
% subplot(6,1,2)
% [h temp] = hist(Time_Peak_Infections,50);
% b = bar(temp,h/nsim);
% set(b,'FaceColor',[51 51 255]/255,'EdgeColor','k','Linewidth',1.5);
% set(gca,'Linewidth',2)
% set(gca,'Linewidth',2)
% axis tight
% str = sprintf('R_0 = %.2f, I2_0=%i (Mean = %.2f, Std = %f)',R0, I2_0, mean(Time_Peak_Infections),std(Time_Peak_Infections));
% title(str);
% ylabel('Frequency')
% xlabel('Time to peak infection [days]')
% grid on
% grid minor
% grid minor
% 
% subplot(6,1,3)
% [h temp] = hist(Time_Outbreak,50);
% b = bar(temp,h/nsim);
% set(b,'FaceColor',[51 51 255]/255,'EdgeColor','k','Linewidth',1.5);
% set(gca,'Linewidth',2)
% set(gca,'Linewidth',2)
% axis tight
% str = sprintf('R_0 = %.2f, I2_0=%i (Mean = %.2f, Std = %f)',R0, I2_0, mean(Time_Outbreak),std(Time_Outbreak));
% title(str);
% ylabel('Frequency')
% xlabel('Time to outbreak [days]')
% grid on
% grid minor
% grid minor
% 
% subplot(6,1,4)
% [h temp] = hist(Num_Deaths,50);
% b = bar(temp,h/nsim);
% set(b,'FaceColor',[51 51 255]/255,'EdgeColor','k','Linewidth',1.5);
% set(gca,'Linewidth',2)
% set(gca,'Linewidth',2)
% axis tight
% str = sprintf('R_0 = %.2f, I2_0=%i (Mean = %.2f, Std = %f)',R0, I2_0, mean(Num_Deaths),std(Num_Deaths));
% title(str);
% ylabel('Frequency')
% xlabel('Number of deaths')
% grid on
% grid minor
% grid minor
% 
% subplot(6,1,5)
% [h temp] = hist(Num_Inf,50);
% b = bar(temp,h/nsim);
% set(b,'FaceColor',[51 51 255]/255,'EdgeColor','k','Linewidth',1.5);
% set(gca,'Linewidth',2)
% set(gca,'Linewidth',2)
% axis tight
% str = sprintf('R_0 = %.2f, I2_0=%i (Mean = %.2f, Std = %f)',R0, I2_0, mean(Num_Inf),std(Num_Inf));
% title(str);
% ylabel('Frequency')
% xlabel('Number of infections')
% grid on
% grid minor
% grid minor
% 
% subplot(6,1,6)
% [h temp] = hist(Num_SilInf,50);
% b = bar(temp,h/nsim);
% set(b,'FaceColor',[51 51 255]/255,'EdgeColor','k','Linewidth',1.5);
% set(gca,'Linewidth',2)
% set(gca,'Linewidth',2)
% axis tight
% str = sprintf('R_0 = %.2f, I2_0=%i (Mean = %.2f, Std = %f)',R0, I2_0, mean(Num_SilInf),std(Num_SilInf));
% title(str);
% ylabel('Frequency')
% xlabel('Number of infected silent spreaders')
% grid on
% grid minor
% grid minor
% 
% str = sprintf(['Hybrid_COVID19_NumDeaths_PeakInf_TimePeak_NumInf_TimeOutbreak_NumSilentInf_r',r_value_print,'_CV',CV_value_print,'_InitialCondition',num2str(InitialCondition)]);
% set(gcf,'Name',str,'NumberTitle','off')
% saveas(gcf,strcat(str,'.fig'))
% saveas(gcf,strcat(str,'.png'))

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
disp((sprintf('R0 = %.2f',R0)))
disp((sprintf('Probability of outbreak = %.2f',numOutbreak/nsim)))
disp((sprintf('Probability of ICU overload = %.2f',numICUoverload/nsim)))
disp((sprintf('----------------------------------------------------------------------')))
disp((sprintf('PeakInf \t %.2f \t %.2f',mean(Peak_Infections),std(Peak_Infections))))
disp((sprintf('TimePeakInf \t %.2f \t %.2f',mean(Time_Peak_Infections),std(Time_Peak_Infections))))
disp((sprintf('NumDeaths \t %.2f \t %.2f',mean(Num_Deaths),std(Num_Deaths))))
disp((sprintf('NumInf \t\t %.2f \t %.2f',mean(Num_Inf),std(Num_Inf))))
disp((sprintf('TimeOutbreak \t %.2f \t %.2f',mean(Time_Outbreak),std(Time_Outbreak))))
disp((sprintf('----------------------------------------------------------------------')))

