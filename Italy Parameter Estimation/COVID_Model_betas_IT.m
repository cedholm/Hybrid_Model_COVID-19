%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COVID_Model.m
%Christina Edholm 
%Edited by Karen Hwang
%g2 = 1/14 instead of 1/28
%Fit d1, assume d2 remains as 1/2.3, d1 < d2
%
% This code is the ODE model we are using, first we establish parameters
% then write out the system of ODEs. The SilentSpreaders are the 1.
%
% Here it is assumed that beta_1 is the transmission rate from A_1, beta_2 from A_2 and beta_3 from I_2
%
% June 1, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dydt = f(t,y,z)

%Parameters

%Parameter Values - COVID

b1=z(1);            %transmission rate for SilentSpreader Asymptomatics
b2=z(2);            %transmission rate for SymptomaticSpreader Asymptomatics
b3=z(3);              %transmission rate for SymptomaticSpreader Infectious
aa1=1/(5.5-2.3);    %rate of transition from exposed to asymptomatic stage for SilentSpreader - COVID -- double check
aa2=1/(5.5-2.3);    %rate of transition from exposed to asymptomatic stage for SymptomaticSpreader - COVID -- double check
d1= z(6);          %rate of transition from asymptomatic to infected stage for SilentSpreader - COVID
d2= 1/2.3;          %rate of transition from asymptomatic to infected stage for SymptomaticSpreader - COVID
mI2=z(4);           %disease-induced mortality rate for infected SymptomaticSpreader - MERS -- region dependent
g2=z(7);            %removal rate for SymptomaticSpreader - COVID

%R0=(b1*N10)/(d1*N)+(b2*N20)/(d2*N)+(b3*N20)/((g2+mI2)*N);
%%p= z(4);

%b3= ((R0-(b1*p)/(d1)-(b2*(1-p))/(d2)))*((g2+mI2)))/(1-p);

%b3=((R0-(b1*p)/d1-(b2*(1-p))/d2)*(g2+mI2))/(1-p);

%b3=((R0-((b1*p)/d1)-((b2*(1-p))/d2))*((g2+mI2))/(1-p));



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
cumcases=d2*a2;
deaths=mI2*y2;

dydt=[s1dot; e1dot; a1dot; r1dot; s2dot; e2dot; a2dot; y2dot; r2dot; cumcases; deaths];

%dydt=dydt'

end


