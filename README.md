MULTIVARIATE TESTS FOR TREND DETECTION

The document is a technical guide for a code designed to calculate and compare the performance of four different tests used for detecting trends in multivariate data series. These tests include (MOT/MDT) and (CET/CIT). It explains how data is simulated using Kendall's Tau coefficients and a copula function, followed by p-value calculations for each test. The document includes instructions on setting up the code, choosing methods for comparison, and details on the inputs and outputs of the program. It specifically advises on installing required packages, setting up an Excel file for results, and handling potential errors or issues during the process.
These codes calculate the performance of the following four tests:
1.	Multivariate overall trend (MOT) test
2.	Multivariate dependence trend (MDT) test
3.	Multivariate Covariance Inversion Test (CIT) 
4.	Multivariate Covariance Eigen-value Test (CET)
Instructions:
-	Please download the GitHub repository containing the code into a local folder and run each script (P-value, Statistique_tests, Power evaluation_trend_Dependence, Power evaluation_ trend tests_Margins, Power evaluation_tests_Dependence_Margins) in RStudio to make the functions available. Ensure that you have the specified packages installe("copula","Kendall","VGAM","gtools","openxlsx","resample","openxlsx","foreach","MASS")                                                        
-	First: run the excel_setup() function once with the directory of your choice.
-	Then, you can launch the different function with the desired inputs. The run_num input determines the placement of the simulation results in the Excel table.
-	The first simulation should be run with run_num = 1, the second with run_num = 2, etc.
-	If the number is not advanced, the result will simply replace the data on the line.
Run excel_setup() ONLY ONCE. If you run the function again with the same input directory, it will display an error. To resolve this, remove the "COMPARE.xlsx" file.
NEVER have the Excel document "COMPARE.xlsx" open when the different functions (is running. This will cause an error at the end of the simulation and the results will be lost.
The coefs input should have more than 20 elements: a non-critical error may occur during calculation.
The inputs:
run_num: A number indicating the placement of results in the Excel table.
•	If run_num = 0, the function will not record values in the Excel file.
•	copula: A string that specifies the copula family. The following families are implemented:
•	Archimedean: clayton, gumbel, joe, frank
•	EV: galambos, huslerReiss
•	coefs: A vector (m x 1) that provides the association coefficients. Values must be compatible with the chosen copula family. A sequence of the type:
•	seq(0,1,length.out=m): Data series with trends in dependency (Power evaluation: trend in dependence)
•	rep(tau,m): Data series without trends in dependency (First type error estimates)
•	m: sample size (30 or 50 or 100)
•	Nsim: Number of simulations (default = 1000)
•	Nbs: Number of resampling for bootstrap estimation of p-value (default = 1000)
•	width: Width of rolling window for the MOT/MDT Tests (default = 10),(10 when m=30, 15 when m=50 and 20 when m=100)
•	alpha: Significance level of the test (default = 0.05)
•	Location_param_x : trend in the first marginal distribution X
•	Location_param_y : trend in the second marginal distribution Y 
•	Trend_in_margin = "both" or "margin_1"

The outputs 
RATES: A vector of length 4 containing the percentages of p-values less than alpha for each method. If a method is not used, a value of 0 is returned.
The function creates an Excel file "COMPARE.xlsx" with the results of previous runs and those of the current run.
     Examples: See Exemples.RMD file
