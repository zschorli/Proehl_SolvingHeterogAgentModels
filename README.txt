
This folder contains MATLAB software accompanying the working
paper “Approximating Equilibria with Ex-Post Heterogeneity
and Aggregate Risk” by Elisabeth Proehl, (October 2018). 

This version: October 2018
The following items are provided: 

1. LICENSE AGREEMENT.

2. FOLDERS. 

   a. Folders containing novel code written for this working paper:

	(1) "Proehl_GrowthModel": Implements the proximal point algorithm 
            and the policy function iteration using polynomial chaos 
            expansions to solve the Aiyagari-Bewley economy, 
            computes stationary distributions and errors, produces plots
	

	(2) "Proehl_HuggettModel": Implements the proximal point algorithm 
            and the policy function iteration using polynomial chaos 
            expansions to solve the Hugget economy, 
            computes stationary distributions and errors, produces plots

   b. Folders containing code from existing papers for the comparison in 
      this working paper: 

      The code in all of the folowing folders was downloaded from 
      http://www.wouterdenhaan.com/datasuite.htm. The code was then
      modified to set the appropriate grids and model parameters for the
      code comparison.

	(3) "DenHaanRendahl": Contains code from the paper
            Wouter J. Den Haan, Pontus Rendahl, Solving the incomplete 
            markets model with aggregate uncertainty using explicit 
            aggregation, Journal of Economic Dynamics and Control, 
            Volume 34, Issue 1, 2010, Pages 69-78.		

	(4) "KrusellSmithByMaliarMaliarValli_1mom": Contains code from 
            the paper
            Lilia Maliar, Serguei Maliar, Fernando Valli, Solving the 
            incomplete markets model with aggregate uncertainty using 
            the Krusell–Smith algorithm, Journal of Economic Dynamics and 
            Control, Volume 34, Issue 1, 2010, Pages 42-49.

	(5) "KrusellSmithByMaliarMaliarValli_2mom": Same as (2) but 
            modified to accommodate second moments

	(6) "KrusellSmithByMaliarMaliarValli_3mom": Same as (2) but 
            modified to accommodate second and third moments

	(7) "KrusellSmithByMaliarMaliarValli_4mom": Same as (2) but 
            modified to accommodate second, third and fourth moments

	(8) "Reiter": Contains code from the paper
	    Michael Reiter, Solving the incomplete markets model with 
            aggregate uncertainty by backward induction, Journal of 
            Economic Dynamics and Control, Volume 34, Issue 1, 2010, 
            Pages 28-35.

3. INSTRUCTIONS TO REPRODUCE ALL RESULTS IN THE PAPER.
	
   a. Growth model:
	
	(1) In folder (1): Set case_nr equal to one in line 36 of "main.m" and 
	    comment out the results section in "main.m". Then run "main.m".

	(2) In folder (3): Define the appropriate path to folder (1) in 
	    line 105 of "ExplicitAggr.m". Run "ExplicitAggr.m" to produce 
	    the result files "DR_Sol1.mat" to "DR_Sol4.mat". Copy the result 
	    file ending with number "Sol..." to 
	    "Proehl_GrowthModel/res_case...".

	(3) In folder (4): Define the appropriate path to folder (1) in 
	    line 133 of "MAIN.m". Run "MAIN.m" to produce the result
	    files "KS1_Sol1.mat" to "KS1_Sol4.mat". Copy the result 
	    file ending with number "Sol..." to 
	    "Proehl_GrowthModel/res_case...".

	(4) In folder (5): Define the appropriate path to folder (1) in 
	    line 127 of "MAIN.m". Run "MAIN.m" to produce the result
	    files "KS2_Sol1.mat" to "KS2_Sol4.mat". Copy the result 
	    file ending with number "Sol..." to 
	    "Proehl_GrowthModel/res_case...".

	(5) In folder (6): Define the appropriate path to folder (1) in 
	    line 127 of "MAIN.m". Run "MAIN.m" to produce the result
	    files "KS3_Sol1.mat" to "KS3_Sol4.mat". Copy the result 
	    file ending with number "Sol..." to 
	    "Proehl_GrowthModel/res_case...".

	(6) In folder (7): Define the appropriate path to folder (1) in 
	    line 127 of "MAIN.m". Run "MAIN.m" to produce the result
	    files "KS4_Sol1.mat" to "KS4_Sol4.mat". Copy the result 
	    file ending with number "Sol..." to 
	    "Proehl_GrowthModel/res_case...".

	(7) In folder (8): Define the appropriate path to folder (1) in 
	    line 191 of "setparam.m". Run "main.m" to produce the result
	    files "R_Sol1.mat" to "R_Sol4.mat". Copy the result 
	    file ending with number "Sol..." to 
	    "Proehl_GrowthModel/res_case...".

	(8) In folder (1): Run the results section in "main.m".

	(9) In folder (1): Vary case_nr in line 36 of "main.m" to run the 
	    different model configurations for the robustness checks. 

	(10) In folder (1): Run "plots.m".

   b. Huggett model:

	(1) In folder (2): Set case_nr equal to one in line 29 of "main.m" 
	    and run "main.m".

	(2) In folder (2): Vary case_nr in line 29 of "main.m" to run the 
	    different model configurations for the robustness checks. 

	(3) In folder (2): Run "plots.m".
   


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
Copyright © 2018 by Elisabeth Proehl. All rights reserved. The code may 
be used, modified and redistributed under the terms provided in the file 
"LICENSE.txt".
-------------------------------------------------------------------------
