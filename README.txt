
This folder contains MATLAB software accompanying the working
paper “Approximating Equilibria with Ex-Post Heterogeneity
and Aggregate Risk” by Elisabeth Proehl, (June 2025). 

This version: June 2025
The following items are provided: 

1. LICENSE AGREEMENT.

2. VIDEOS OF SIMULATIONS

3. CODE CONTAINING SUBFOLDERS:

   a. Folders containing novel code written for this working paper:

	(1) "Proehl_GrowthModel": Implements the proximal point algorithm 
            and the policy function iteration using polynomial chaos 
            expansions to solve the Aiyagari-Bewley economy with aggregate  
            risk, computes stationary distributions and errors, produces plots
	

	(2) "Proehl_EndLaborGrowthModel": Implements the policy function 
			iteration using polynomial chaos expansions to solve
            the Aiyagari-Bewley economy with aggregate risk and  
            endogenous labor supply choice, computes stationary 
			distributions and errors, produces plots/tables

   b. Folders containing code from existing papers for the comparison in 
      this working paper: 

      The code in all of the folowing folders was downloaded from 
      http://www.wouterdenhaan.com/datasuite.htm. The code was then
      modified to set the appropriate grids and model parameters for the
      code comparison.

	(3) "KrusellSmithByMaliarMaliarValli_1mom": Contains code from 
            the paper
            Lilia Maliar, Serguei Maliar, Fernando Valli, Solving the 
            incomplete markets model with aggregate uncertainty using 
            the Krusell–Smith algorithm, Journal of Economic Dynamics and 
            Control, Volume 34, Issue 1, 2010, Pages 42-49.

	(4) "KrusellSmithByMaliarMaliarValli_2mom": Same as (2) but 
            modified to accommodate second moments

	(5) "KrusellSmithByMaliarMaliarValli_3mom": Same as (2) but 
            modified to accommodate second and third moments

	(6) "KrusellSmithByMaliarMaliarValli_4mom": Same as (2) but 
            modified to accommodate second, third and fourth moments


4. INSTRUCTIONS TO REPRODUCE ALL RESULTS IN THE PAPER.
	
   a. Baseline model:
	
	1. 	In Folders (3)-(6): Define the appropriate path to Folder 1 in lines 
		116 and 122 of ”MAIN.m”.

	2. 	In Folder (1): Define the appropriate path to Folders 3-6 in line 61
		of ”main.m”.

	3. 	In Folder (1): Define the appropriate path to Folder 1 in line 70 of
		”main.m”.

	4. 	In Folder (1): Run ”main.m”.

	5. 	In Folder (1): Run ”plots.m” to reproduce the figures.

   b. Extended model:

	1. 	In Folder (2): Run ”main.m”.

	2.	In Folder (2): Run ”plots.m” to reproduce the figures and tables.
   

For updates, please check the authors' web page 
(www.elisbethproehl.com). For additional information, please 
contact the author: 

Elisabeth Proehl
University of Amsterdam
Roeterstraat 11
1018 WB Amsterdam
Netherlands

e.proehl@uva.nl

-------------------------------------------------------------------------
Copyright © 2025by Elisabeth Proehl. All rights reserved. The code may 
be used, modified and redistributed under the terms provided in the file 
"LICENSE.txt".
-------------------------------------------------------------------------
