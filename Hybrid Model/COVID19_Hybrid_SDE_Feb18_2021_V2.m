
function [timestep,time_curr,beta1_overtime,beta2_overtime,beta3_overtime,S1_overtime,S2_overtime,E1_overtime,E2_overtime,A1_overtime,A2_overtime,I2_overtime,R1_overtime,R2_overtime,numDeaths,Switch] = COVID19_Hybrid_SDE(timestep,time_curr,beta1_overtime,beta2_overtime,beta3_overtime,S1_overtime,S2_overtime,E1_overtime,E2_overtime,A1_overtime,A2_overtime,I2_overtime,R1_overtime,R2_overtime,numDeaths,threshold_hybrid_SDEtoNHP,Switch);
   
    global fraction_silent fraction_symp R0_fit b1_0 b2_0 b3_0 mI2 aa1 aa2 dd1 dd2 g2 POPULATION N_0 N1_0 N2_0 R0 r1 r2 r3 sigma1 sigma2 sigma3 time dt
    
thresh=.01; %Don't change (ensure solutions nonnegative when variables small)

%% Time parameters & number of simulatioms
sdt=sqrt(dt);
ts=time/dt;
nsim=1;
%% Beginning 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial Conditions for COVID Model 
%
R1=R1_overtime(timestep);              %Susceptible individuals for SilentSpreader
E1=E1_overtime(timestep);              %Exposed individuals but not infectious for SilentSpreader
A1=A1_overtime(timestep);              %Asymptomatic individuals but infectious for SilentSpreader
%
R2=R2_overtime(timestep);              %Susceptible individuals for SymptomaticSpreaders
E2=E2_overtime(timestep);              %Exposed individuals but not infectious for SymptomaticSpreaders
I2=I2_overtime(timestep);              %Infected individuals for SymptomaticSpreaders
A2=A2_overtime(timestep);              %Asymptomatic individuals but infectious for SymptomaticSpreader
%
S1=S1_overtime(timestep);            %SilentSpreaders
S2=S2_overtime(timestep);        %SymptomaticSpreaders


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%Euler Maruyama  method for solving Full SDE 
%Define initial conditions for numerical  method as VECTOR
%
         beta1=beta1_overtime(timestep);
         beta2=beta2_overtime(timestep);
         beta3=beta3_overtime(timestep);        
         Deaths=zeros(nsim,1);   %NEW
         
while E1+E2+A1+A2+I2>0 && time_curr<time && Switch == 1
        rnn=randn(nsim,8)*sdt; 
        %Total Population
        N10=S1+E1+A1+R1;
        N20=S2+E2+A2+I2+R2;       
%Make all rates nonnegative (do not change the values of random 
%variables S, E, A. I, or R, just the rates!) terms with 11 or 22  
% ensure nonnegative rates same as the CTMC model.
% Terms with 'w' are positive, as they are in the Wiener process
        fS11=(S1./(N10+N20)).*(beta1.*A1+beta2.*A2+beta3.*I2); 
        fS11pos=fS11>=0; fS11=fS11pos.*fS11;
        fS1=-fS11;
        fS1w=fS11; 
        fE1a= fS11;
        %
        fE11b=aa1*E1;
        fE11bpos=fE11b>=0; fE11b=fE11bpos.*fE11b;
        fE1b=-fE11b;
        fE1bw=fE11b; 
        fA1a=fE11b;
        %
        fA11b=dd1*A1;
        fA11bpos=fA11b>=0; fA11b=fA11bpos.*fA11b;
        fA1b=-fA11b;
        fA1bw=fA11b; 
        fR1=fA11b;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fS22=(S2./(N10+N20)).*(beta1.*A1+beta2.*A2+beta3.*I2);
        fS22pos=fS22>=0; fS22=fS22pos.*fS22;
        fS2=-fS22;
        fS2w=fS22; 
        fE2a= fS22;
        %
        fE22b=aa2*E2;
        fE22bpos=fE22b>=0; fE22b=fE22bpos.*fE22b;
        fE2b=-fE22b;
        fE2bw=fE22b; 
        fA2a=fE22b;
        %
        fA22b=dd2*A2;
        fA22bpos=fA22b>=0; fA22b=fA22bpos.*fA22b;
        fA2b=-fA22b;
        fA2bw=fA22b; 
        fI2a=fA22b; 

        fI22b=mI2*I2;
        fI22bpos=fI22b>=0; fI22b=fI22bpos.*fI22b;
        fI2b=-fI22b;
        fI2bw=fI22b;
        %
        fI22c=g2*I2;
        fI22cpos=fI22c>=0; fI22c=fI22cpos.*fI22c;
        fI2c=-fI22c;
        fI2cw=fI22c;  
        fR2=fI22c;
        %
        %Thresh is the cutoff value where
        %diffusion terms are modified below thresh.  This
        %step ensures the SDE has nonnegative solutions. 
        %Only changes 3 WIENER Process rates in silent and 4 in symptomatic (total 7).
        %
        % Makes diffusion terms continuous 
          E1thresh=(E1<=thresh);E1pos=(E1>=0);
          fS1w=E1thresh.*fS1w.*E1pos.*E1/thresh+(ones(nsim,1)-E1thresh).*fS1w;
          A1thresh=(A1<=thresh);A1pos=(A1>=0);
          fE1bw=A1thresh.*fE1bw.*A1pos.*A1/thresh+(ones(nsim,1)-A1thresh).*fE1bw;
          %
         % A1thresh=(A1<=thresh);A1pos=(A1>=0); %%not needed for COVID
         % fA1bw=A1thresh.*fA1bw.*A1pos.*A1/thresh+(ones(nsim,1)-A1thresh).*fA1bw;
         %
         %I1thresh=(I1<=thresh);I1pos=(I1>=0);
         %fA1bw=I1thresh.*fA1bw.*I1pos.*I1/thresh+(ones(nsim,1)-I1thresh).*fA1bw;
          R1thresh=(R1<=thresh); R1pos=(R1>=0);
          fA1bw=R1thresh.*fA1bw.*R1pos.*R1/thresh+(ones(nsim,1)-R1thresh).*fA1bw;
          %
          E2thresh=(E2<=thresh); E2pos=(E2>=0);
          fS2w=E2thresh.*fS2w.*E2pos.*E2/thresh+(ones(nsim,1)-E2thresh).*fS2w;
          A2thresh=(A2<=thresh);A2pos=(A2>=0);
          fE2bw=A2thresh.*fE2bw.*A2pos.*A2/thresh+(ones(nsim,1)-A2thresh).*fE2bw;
          I2thresh=(I2<=thresh); I2pos=(I2>=0);
          fA2bw=I2thresh.*fA2bw.*I2pos.*I2/thresh+(ones(nsim,1)-I2thresh).*fA2bw;
          R2thresh=(R2<=thresh);R2pos=(R2>=0);
          fI2cw=R2thresh.*fI2cw.*R2pos.*R2/thresh+(ones(nsim,1)-R2thresh).*fI2cw;
           %
       % Solve the SDE system with Euler-Maruyama
       % Note here we do not need absolute value 
       % because all rates are positive and some 'w' terms are
       % revised as above so that solutions are nonnegative.
       %
        S1= S1+dt*fS1-sqrt(fS1w).*rnn(:,1);
        S2= S2+dt*fS2-sqrt(fS2w).*rnn(:,4);
        %
        E1= E1+dt*(fE1a+fE1b)+sqrt(fS1w).*rnn(:,1)-sqrt(fE1bw).*rnn(:,2);
        E2= E2+dt*(fE2a+fE2b)+sqrt(fS2w).*rnn(:,4)-sqrt(fE2bw).*rnn(:,5);
%
        A1= A1+dt*(fA1a+fA1b)+sqrt(fE1bw).*rnn(:,2)-sqrt(fA1bw).*rnn(:,3); % -sqrt(fA1cw).*rnn(:,4);
        A2= A2+dt*(fA2a+fA2b)+sqrt(fE2bw).*rnn(:,5)-sqrt(fA2bw).*rnn(:,6); %-sqrt(fA2cw).*rnn(:,10);
%
       % I1= I1+dt*(fI1a+fI1b+fI1c)+sqrt(fA1bw).*rnn(:,3)-sqrt(fI1bw).*rnn(:,5)...
       %        -sqrt(fI1cw).*rnn(:,6); %No I1 compartment in the COVID model 
       
        I2= I2+dt*(fI2a+fI2b+fI2c)+sqrt(fA2bw).*rnn(:,6)-sqrt(fI2bw).*rnn(:,7)...
                -sqrt(fI2cw).*rnn(:,8);
            
        Deaths=Deaths-fI2b*dt+sqrt(fI2bw).*rnn(:,7); %%NEW
        
        R1= R1+dt*(fR1)+sqrt(fA1bw).*rnn(:,3);
        R2= R2+dt*(fR2)+sqrt(fI2cw).*rnn(:,8);
        %
        % SDE for environmental variability
        rn=randn(nsim,1)*sdt;
        %
        beta1=max(0,beta1+dt*r1*(b1_0-beta1)+sigma1*rn(:,1).*sqrt(beta1));
        beta2=max(0,beta2+dt*r2*(b2_0-beta2)+sigma2*rn(:,1).*sqrt(beta2));   
        beta3=max(0,beta3+dt*r3*(b3_0-beta3)+sigma3*rn(:,1).*sqrt(beta3));
        %

        if E1+E2+A1+A2+I2<=threshold_hybrid_SDEtoNHP
            Switch = 0;
        end
        
        time_curr=time_curr+dt;
        timestep=timestep+1;
        
        S1_overtime(timestep)=S1;
        S2_overtime(timestep)=S2;
        E1_overtime(timestep)=E1;
        E2_overtime(timestep)=E2;
        A1_overtime(timestep)=A1;
        A2_overtime(timestep)=A2;
        I2_overtime(timestep)=I2;
        R1_overtime(timestep)=R1;
        R2_overtime(timestep)=R2;
        beta1_overtime(timestep)=beta1;
        beta2_overtime(timestep)=beta2;
        beta3_overtime(timestep)=beta3;
        
    end
    
    numDeaths=Deaths+numDeaths;
end