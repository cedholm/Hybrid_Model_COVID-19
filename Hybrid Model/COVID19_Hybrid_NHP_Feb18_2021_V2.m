
% Written by: Karen K. L. Hwang, Christina J. Edholm, Omar Saucedo, Linda J. S. Allen, Nika Shakiba

function [timestep,time_curr,beta1_overtime,beta2_overtime,beta3_overtime,S1_overtime,S2_overtime,E1_overtime,E2_overtime,...
    A1_overtime,A2_overtime,I2_overtime,R1_overtime,R2_overtime,numDeaths,Error_check,Switch] =...
    COVID19_Hybrid_NHP(timestep,time_curr,kk,beta1_overtime,beta2_overtime,beta3_overtime,S1_overtime,S2_overtime,E1_overtime,...
    E2_overtime,A1_overtime,A2_overtime,I2_overtime,R1_overtime,R2_overtime,numDeaths,Error_check,threshold_hybrid_NHPtoSDE,Switch)

    global fraction_silent fraction_symp R0_fit b1_0 b2_0 b3_0 mI2 aa1 aa2 dd1 dd2 g2 POPULATION N_0 N1_0 N2_0 R0 r1 r2 r3 sigma1 sigma2 sigma3 time dt

    sdt = sqrt(dt);

    t = 1;
    
    beta1(t) = beta1_overtime(timestep);
    beta2(t) = beta2_overtime(timestep);
    beta3(t) = beta3_overtime(timestep);
    E1(t)=E1_overtime(timestep);            E2(t)=E2_overtime(timestep);
    A1(t)=A1_overtime(timestep);            A2(t)=A2_overtime(timestep);
                                            I2(t)=I2_overtime(timestep);
    R1(t)=R1_overtime(timestep);            R2(t)=R2_overtime(timestep);
    S1(t)=S1_overtime(timestep);            S2(t)=S2_overtime(timestep);
    
    while E1(t)+E2(t)+A1(t)+A2(t)+I2(t)>0 && time_curr<time && Switch == 0
        
        trans=(beta1(t))*(A1(t))+(beta2(t))*(A2(t))+(beta3(t))*(I2(t));
        totN=S1(t)+S2(t)+E1(t)+E2(t)+A1(t)+A2(t)+I2(t)+R1(t)+R2(t);
        ev1=S1(t)*trans/totN*dt;                                            % S1 --> E1
        ev2=ev1+S2(t)*trans/totN*dt;                                        % S2 --> E2
        ev3=ev2+aa1*E1(t)*dt;                                               % E1 --> A1
        ev4=ev3+aa2*E2(t)*dt;                                               % E2 --> A2
        ev5=ev4+dd1*A1(t)*dt;                                               % A1 --> R1
        ev6=ev5+dd2*A2(t)*dt;                                               % A2 --> I2
        ev7=ev6+mI2*I2(t)*dt;                                               % I2 --> death
        ev8=ev7+g2*I2(t)*dt;                                                % I2 --> R2

        %Check whether the probabilities<1
        %if ev8>1 or summ>0 the time step dt needs to be smaller
        if ev8>1
            Error_check(kk,timestep) = Error_check(kk,timestep) + 1;
        end

        rn=rand; %uniform random number generator
        rann=randn; %normal random number generator
        rnn=rann*sdt;
        
        if rn<ev1
            S1(t+1)=S1(t)-1;
            S2(t+1)=S2(t);
            E1(t+1)=E1(t)+1;
            E2(t+1)=E2(t);
            A1(t+1)=A1(t);
            A2(t+1)=A2(t);
            I2(t+1)=I2(t);
            R1(t+1)=R1(t);
            R2(t+1)=R2(t);
        elseif rn>=ev1 && rn<ev2
            S1(t+1)=S1(t);
            S2(t+1)=S2(t)-1;
            E1(t+1)=E1(t);
            E2(t+1)=E2(t)+1;
            A1(t+1)=A1(t);
            A2(t+1)=A2(t);
            I2(t+1)=I2(t);
            R1(t+1)=R1(t);
            R2(t+1)=R2(t);
        elseif rn>=ev2 && rn<ev3
            S1(t+1)=S1(t);
            S2(t+1)=S2(t);
            E1(t+1)=E1(t)-1;
            E2(t+1)=E2(t);
            A1(t+1)=A1(t)+1;
            A2(t+1)=A2(t);
            I2(t+1)=I2(t);
            R1(t+1)=R1(t);
            R2(t+1)=R2(t);
        elseif rn>=ev3 && rn<ev4
            S1(t+1)=S1(t);
            S2(t+1)=S2(t);
            E1(t+1)=E1(t);
            E2(t+1)=E2(t)-1;
            A1(t+1)=A1(t);
            A2(t+1)=A2(t)+1;
            I2(t+1)=I2(t);
            R1(t+1)=R1(t);
            R2(t+1)=R2(t);
        elseif rn>=ev4 && rn<ev5
            S1(t+1)=S1(t);
            S2(t+1)=S2(t);
            E1(t+1)=E1(t);
            E2(t+1)=E2(t);
            A1(t+1)=A1(t)-1;
            A2(t+1)=A2(t);
            I2(t+1)=I2(t);
            R1(t+1)=R1(t)+1;
            R2(t+1)=R2(t);
        elseif rn>=ev5 && rn<ev6
            S1(t+1)=S1(t);
            S2(t+1)=S2(t);
            E1(t+1)=E1(t);
            E2(t+1)=E2(t);
            A1(t+1)=A1(t);
            A2(t+1)=A2(t)-1;
            I2(t+1)=I2(t)+1;
            R1(t+1)=R1(t);
            R2(t+1)=R2(t);
        elseif rn>=ev6 && rn<ev7
            S1(t+1)=S1(t);
            S2(t+1)=S2(t);
            E1(t+1)=E1(t);
            E2(t+1)=E2(t);
            A1(t+1)=A1(t);
            A2(t+1)=A2(t);
            I2(t+1)=I2(t)-1;
            numDeaths = numDeaths + 1;
            R1(t+1)=R1(t);
            R2(t+1)=R2(t);
        elseif rn>=ev7 && rn<ev8
            S1(t+1)=S1(t);
            S2(t+1)=S2(t);
            E1(t+1)=E1(t);
            E2(t+1)=E2(t);
            A1(t+1)=A1(t);
            A2(t+1)=A2(t);
            I2(t+1)=I2(t)-1;
            R1(t+1)=R1(t);
            R2(t+1)=R2(t)+1;
        else
            S1(t+1)=S1(t);
            S2(t+1)=S2(t);
            E1(t+1)=E1(t);
            E2(t+1)=E2(t);
            A1(t+1)=A1(t);
            A2(t+1)=A2(t);
            I2(t+1)=I2(t);
            R1(t+1)=R1(t);
            R2(t+1)=R2(t);
        end

        %to ensure values do not get negative we bound below by zero by
        %taking maximum value 
        %constant beta let constbeta=0
        beta1(t+1)=max(0,beta1(t)+r1*(b1_0-beta1(t))*dt+sigma1*sqrt(abs(beta1(t)))*rnn);
        beta2(t+1)=max(0,beta2(t)+r2*(b2_0-beta2(t))*dt+sigma2*sqrt(abs(beta2(t)))*rnn);
        beta3(t+1)=max(0,beta3(t)+r3*(b3_0-beta3(t))*dt+sigma3*sqrt(abs(beta3(t)))*rnn);

        S1(t)=S1(t+1);
        S2(t)=S2(t+1); 
        E1(t)=E1(t+1);
        E2(t)=E2(t+1);
        A1(t)=A1(t+1);
        A2(t)=A2(t+1);
        I2(t)=I2(t+1);
        R1(t)=R1(t+1);
        R2(t)=R2(t+1);
        beta1(t)=beta1(t+1);
        beta2(t)=beta2(t+1);
        beta3(t)=beta3(t+1);

        timestep = timestep + 1;
        time_curr = dt+time_curr;
        
        S1_overtime(timestep)=S1(t+1);
        S2_overtime(timestep)=S2(t+1); 
        E1_overtime(timestep)=E1(t+1);
        E2_overtime(timestep)=E2(t+1);
        A1_overtime(timestep)=A1(t+1);
        A2_overtime(timestep)=A2(t+1);
        I2_overtime(timestep)=I2(t+1);
        R1_overtime(timestep)=R1(t+1);
        R2_overtime(timestep)=R2(t+1);
        beta1_overtime(timestep)=beta1(t+1);
        beta2_overtime(timestep)=beta2(t+1);
        beta3_overtime(timestep)=beta3(t+1);
                
        if E1(t)+E2(t)+A1(t)+A2(t)+I2(t)>=threshold_hybrid_NHPtoSDE
            Switch = 1;
        end
        
    end

end
