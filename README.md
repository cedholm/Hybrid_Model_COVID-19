# Hybrid_Model_COVID-19

## MATLAB code contents
### 1. Italy parameter estimation ### 
<br /> &nbsp;&nbsp;&nbsp;
[Italy Days 1-22](https://github.com/cedholm/Hybrid_Model_COVID-19/tree/main/Italy%20Parameter%20Estimation)
<br /> &nbsp;&nbsp;&nbsp;
Used to generate parameter fitting results for Italy from Days 1-22 (February 20, 2020 - March 12, 2020)
<br /> &nbsp;&nbsp;


### 2. BC parameter estimation ### 
<br /> &nbsp;&nbsp;&nbsp;
[BC Phase I](https://github.com/cedholm/Hybrid_Model_COVID-19/tree/main/BC%20Parameter%20Estimation/Phase%201)
<br /> &nbsp;&nbsp;
Used to generate parameter fitting results for BC from Days 1-60 (January 30, 2020 - March 31, 2020)
<br /> &nbsp;&nbsp;
<br /> &nbsp;&nbsp;
[BC Phase II](https://github.com/cedholm/Hybrid_Model_COVID-19/tree/main/BC%20Parameter%20Estimation/Phase%202)
<br /> &nbsp;&nbsp;
Used to generate parameter fitting results for BC from Days 60-109 (March 31, 2020 - May 9, 2020)
<br /> &nbsp;&nbsp;
<br /> &nbsp;&nbsp;


## How to run MultiStart
The main file name is COVID_MultiStart.m, which calls COVID_model.m. Fixed and fitted parameter values are specific to the version of the code included in each section. The fixed parameters are listed in COVID_model.m.
<br /> &nbsp;&nbsp;

### Inputs ###
<br /> &nbsp;&nbsp;&nbsp;
*COVID_MultiStart(NoStartPoints, Tstart, Tend, place, testnumber)*
<br /> &nbsp;&nbsp;&nbsp;
- **NoStartPoints**: The number of starting points you want MultiStart to use, the more you use the longer the code will take to run. 
- **Tstart**: Day of data you want to start fitting the model to 
- **Tend**: Day of data you want to end fitting the model to 
- **Place**: Location you are fitting for ('BC' for British Columbia or 'IT' for Italy)
- **Testnumber**: Test number recorded on all outputs to help you keep track of then number of runs
<br /> &nbsp;&nbsp;&nbsp;
<br /> &nbsp;&nbsp;&nbsp;
*Note: Ensure the corresponding excel data file is in the same folder as both COVID_MultiStart.m and COVID_model.m.*
<br /> &nbsp;&nbsp;&nbsp;
