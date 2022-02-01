# Hybrid_Model_COVID-19

## MATLAB code contents
### 1. Italy parameter estimation ### 
<br /> &nbsp;&nbsp;
[Italy Days 1-22](https://github.com/cedholm/Hybrid_Model_COVID-19/tree/main/Italy%20Parameter%20Estimation)
<br /> &nbsp;&nbsp;&nbsp;
Used to generate parameter fitting results for Italy from Days 1-22 (February 20, 2020 - March 12, 2020)
<br /> &nbsp;&nbsp;


### 2. BC parameter estimation ### 
<br /> &nbsp;&nbsp;
[BC Phase I](https://github.com/cedholm/Hybrid_Model_COVID-19/tree/main/BC%20Parameter%20Estimation/Phase%201)
<br /> &nbsp;&nbsp;
Used to generate parameter fitting results for BC from Days 1-60 (January 30, 2020 - March 31, 2020)
<br /> &nbsp;&nbsp;
<br /> &nbsp;&nbsp;
[BC Phase II](https://github.com/cedholm/Hybrid_Model_COVID-19/tree/main/BC%20Parameter%20Estimation/Phase%202)
<br /> &nbsp;&nbsp;
Used to generate parameter fitting results for BC from Days 60-109 (March 31, 2020 - May 19, 2020)
<br /> &nbsp;&nbsp;

### 3. Hybrid stochastic model ### 
<br /> &nbsp;&nbsp;
[Hybrid Model](https://github.com/cedholm/Hybrid_Model_COVID-19/tree/main/Hybrid%20Model)
<br /> &nbsp;&nbsp;
The hybrid model calls both the NHP and SDE. Ensure all files in the folder are downloaded
<br /> &nbsp;&nbsp;&nbsp;
<br /> &nbsp;&nbsp;&nbsp;

## How to run MultiStart
The main file name is COVID_MultiStart.m, which calls COVID_model.m. Fixed and fitted parameter values are specific to the version of the code included in each section. The fixed parameters are listed in COVID_model.m.
<br /> &nbsp;&nbsp;

### Inputs ###
<br /> &nbsp;&nbsp;&nbsp;
*COVID_MultiStart(NoStartPoints, Tstart, Tend, place, testnumber)*
<br /> &nbsp;&nbsp;&nbsp;
- **NoStartPoints**: The number of starting points you want MultiStart to use, the more you use the longer the code will take to run 
- **Tstart**: Day of data you want to start fitting the model to 
- **Tend**: Day of data you want to end fitting the model to 
- **Place**: Location you are fitting for ('BC' for British Columbia or 'IT' for Italy)
- **Testnumber**: Test number recorded on all outputs to help you keep track of then number of runs
<br /> &nbsp;&nbsp;&nbsp;
<br /> &nbsp;&nbsp;&nbsp;
*Note: Ensure the corresponding excel data file is in the same folder as both COVID_MultiStart.m and COVID_model.m.*
<br /> &nbsp;&nbsp;&nbsp;

### Outputs ###
- **COVIDParameters**: a matrix of all the possible local minimum parameter values, size will relate to the NoStartPoints
- **fvalues**: all the objective functional minimum values that correspond to the COVIDParameters
- **ExitFlags**: these values give insight into the convergence of fmincon. We want a 1, and will accept a 2. You will get one for each COVIDParameter and fvalues combo.
- **Endpoints**: The state variables values at Tend

<br /> &nbsp;&nbsp;&nbsp;
*When you run the code, it will also save some files in the same folder labelled with the place and the testnumber:*
- An excel file with the COVIDParameters, values, and ExitFlags
- An excel file with the endpoints 
- A figure file of the output trendline plotted on top of data
- A figure of the ODE solution, outlining the number of individuals in S, E, A, I, R, respectively

*The command window will print out:*
- Case Minimization Value
- Death Minimization Value
- Total Minimization Value
- R0
<br /> &nbsp;&nbsp;&nbsp;
<br /> &nbsp;&nbsp;&nbsp;



## How to run hybrid stochastic model
The parameter values for the hybrid model are derived from parameter estimation. Change the values for r and CV to see the change in transmission dynamics associated with environmental variation. Update nsim to change the number of simulations you plan to run. 

