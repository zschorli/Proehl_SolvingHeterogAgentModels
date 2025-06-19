% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This is the main file
%__________________________________________________________________________

%--------------------------------------------------------------------------
% Choose model configuration
%--------------------------------------------------------------------------
% case 1: benchmark model with end. labor (see online appendix on model 
% parameters)
% case 2: benchmark model but linear taxation, i.e. rho=0
% case 3: constant unempl. benefit nu = 0.4
% case 4: constant nu = 0.4 and linear taxation, i.e. rho=0
% case 5: countercycl. nu = 0.135-0.165
% case 6: procycl. nu = 0.135-0.165
% case 7: constant nu = 0
% case 8: constant nu = 0.015
% case 9: countercycl. nu = 0-0.03 
% case 10: procycl. nu = 0-0.03 
% case 11: countercycl. nu = 0-0.015
% case 12: procycl. nu = 0-0.015 
% case 13: countercycl. nu = 0-0.15
% case 14: procycl. nu = 0-0.15 
%
% case 15-24: borrowing constraint k_min = -5
% case 15: constant nu = 0
% case 16: constant nu = 0.015
% case 17: countercycl. nu = 0-0.015 
% case 18: procycl. nu = 0-0.015  
% case 19: constant nu = 0.15
% case 20: countercycl. nu = 0-0.15
% case 21: procycl. nu = 0-0.15
% case 22: constant nu = 0.0075
% case 23: countercycl. nu = 0.01-0.02
% case 24: procycl. nu = 0.01-0.02
%
% case 25-34: borrowing constraint k_min = 5
% case 25: constant nu = 0
% case 26: constant nu = 0.015
% case 27: countercycl. nu = 0-0.015 
% case 28: procycl. nu = 0-0.015  
% case 29: constant nu = 0.15
% case 30: countercycl. nu = 0-0.15
% case 31: procycl. nu = 0-0.15
% case 32: constant nu = 0.0075
% case 33: countercycl. nu = 0.01-0.02
% case 34: procycl. nu = 0.01-0.02
%
% case 35-44: borrowing constraint k_min = -5 (recession), -4.5 (boom)
% case 35: constant nu = 0
% case 36: constant nu = 0.015
% case 37: countercycl. nu = 0-0.015 
% case 38: procycl. nu = 0-0.015  
% case 39: constant nu = 0.15
% case 40: countercycl. nu = 0-0.15
% case 41: procycl. nu = 0-0.15
% case 42: constant nu = 0.0075
% case 43: countercycl. nu = 0.01-0.02
% case 44: procycl. nu = 0.01-0.02

for case_nr = 1:44
folder = strcat('res_case',num2str(case_nr));        
if exist(folder,'dir') ~= 7
    mkdir(folder);
end

%--------------------------------------------------------------------------
% Parameter definition
%--------------------------------------------------------------------------
StaticParams = struct;
StaticParams = setModelParams(StaticParams,case_nr);
save(strcat(folder,'/StaticParams.mat'),'-struct','StaticParams');
end

%--------------------------------------------------------------------------
% Parallel Computing
%--------------------------------------------------------------------------
% Uncomment and modify the following lines if you want to specify a  
% specific number of CPUs to be used for the computation:
% noWorkers = 20;
% parpool(noWorkers);

%--------------------------------------------------------------------------
% Model without aggregate risk
%--------------------------------------------------------------------------
for case_nr = [1:4,7]
folder = strcat('res_case',num2str(case_nr));   
StaticParams = load(strcat(folder,'/StaticParams.mat'));
StaticParams.l_min = 0;
solveModel_NoAggShock(StaticParams,folder);
end

for case_nr = 1:44
folder = strcat('res_case',num2str(case_nr));   
StaticParams = load(strcat(folder,'/StaticParams.mat'));   
%--------------------------------------------------------------------------
% Generate orthogonal polynomials for polynomial chaos approximation
%--------------------------------------------------------------------------
% first set the order of truncation of the polynomial chaos approximation
if case_nr<=15
    order = 3;
else
    switch case_nr
        case {17,18,35,37,38}
            order = 3;
        otherwise
            order =2;
    end
end
setPoly(StaticParams,folder,order,case_nr);

%--------------------------------------------------------------------------
% Model with aggregate risk
%--------------------------------------------------------------------------
solveModel_AggShock(StaticParams,case_nr,folder,order);

%--------------------------------------------------------------------------
% Compute Results
%--------------------------------------------------------------------------
results_statDistr(StaticParams,case_nr,folder,order);
if case_nr == 1
    results_Errors(StaticParams,folder);
end
end