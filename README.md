# CovariateRCT
Code to replicate "Covariate Selection in Randomized Controlled Trials". https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6063750 
Software: R, Numbers (MacOS).

1. The Folder "Main Scripts" contains the scripts to reproduce all the Tables and Figures in the manuscript.

  - "1 - Univariate Case.R": Reproduces Table 1 and Figure 1
  
  - "2 - MSE unbiased estimate.R": Reproduces Table 2
  
  - "3 - MSE Confidence Intervals.R": Reproduces Table 3
  
  - "4 - MSE Multivariate ExactTest.R": Reproduces Table 4 and Figure 2
  
  - "5 - All Multivariate Tests.R": Reproduces Table 5
  
  - "6 - Power.R": Reproduces Table 6 and Figure 3

  - "7 - Selection rules.R" : Reproduces Table 7 and Figure 4
   
2. The Folder "DGP" contains the settings of the data-generating processes used in the simulations.
   
3. The Folder "Package" contains the functions developed to implement the statistical tests presented in the manuscript to compare the MSE of various statistical models.
   
4. The Folder "Function" contains the scripts that define auxiliary functions for each specific script.
   
5. The Folder "Output" will contain the graphs produced by the scripts.

6. The Folder "Graphs" contains the .numbers file to produce the graph of the power analysis. 
