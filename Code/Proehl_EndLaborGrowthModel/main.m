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
% case 15-21: borrowing constraint k_min = 20
% case 15: constant nu = 0
% case 16: constant nu = 0.015
% case 17: countercycl. nu = 0-0.015 
% case 18: procycl. nu = 0-0.015  
% case 19: constant nu = 0.15
% case 20: countercycl. nu = 0-0.15
% case 21: procycl. nu = 0-0.15

% for case_nr = 1:21
% folder = strcat('res_case',num2str(case_nr));        
% if exist(folder,'dir') ~= 7
%     mkdir(folder);
% end
% 
% %--------------------------------------------------------------------------
% % Parameter definition
% %--------------------------------------------------------------------------
% StaticParams = struct;
% StaticParams = setModelParams(StaticParams,case_nr);
% save(strcat(folder,'/StaticParams.mat'),'-struct','StaticParams');
% end

%--------------------------------------------------------------------------
% Parallel Computing
%--------------------------------------------------------------------------
% Uncomment and modify the following lines if you want to specify a  
% specific number of CPUs to be used for the computation:
% noWorkers = 20;
% parpool(noWorkers);

% for case_nr = [1:4,7]
% folder = strcat('res_case',num2str(case_nr));   
% StaticParams = load(strcat(folder,'/StaticParams.mat'));
% %--------------------------------------------------------------------------
% % Model without aggregate risk
% %--------------------------------------------------------------------------
% StaticParams.l_min = 0;
% solveModel_NoAggShock(StaticParams,folder);
% end

for case_nr = [20,21,5,6,13,14]
folder = strcat('res_case',num2str(case_nr));   
StaticParams = load(strcat(folder,'/StaticParams.mat'));   
%--------------------------------------------------------------------------
% Generate orthogonal polynomials for polynomial chaos approximation
%--------------------------------------------------------------------------
% first set the order of truncation of the polynomial chaos approximation
order = 3;
setPoly(StaticParams,folder,order,case_nr);

%--------------------------------------------------------------------------
% Model with aggregate risk
%--------------------------------------------------------------------------
solveModel_AggShock(StaticParams,case_nr,folder,order);

%--------------------------------------------------------------------------
% Compute Results
%--------------------------------------------------------------------------
results_statDistr(StaticParams,case_nr,folder);
if case_nr == 1
    results_Errors(StaticParams,folder);
end
end