% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This is the main file
%__________________________________________________________________________

%--------------------------------------------------------------------------
% Choose model configuration
%--------------------------------------------------------------------------
% Choose a number between 1 and 13:
% 1 benchmark model (see online appendix on model parameters)
% 2 borrowing constraint of -0.5
% 3 borrowing constraint of -1
% 4 borrowing constraint of -2
% 5 borrowing constraint of -2.5
% 6 subjective discount factor of 0.9
% 7 subjective discount factor of 0.999
% 8 risk aversion of 1
% 9 risk aversion of 3
% 10 borrowing constr. -0.5, subjective discount factor of 0.9
% 11 borrowing constr. -0.5, subjective discount factor of 0.999
% 12 borrowing constr. -0.5, risk aversion of 1
% 13 borrowing constr. -0.5, risk aversion of 3

case_nr = 1;

folder = strcat('res_Huggett_case',num2str(case_nr));        
if exist(folder,'dir') ~= 7
    mkdir(folder);
end

%--------------------------------------------------------------------------
% Parameter definition
%--------------------------------------------------------------------------
StaticParams = struct;
StaticParams = setModelParams(StaticParams,case_nr);
save(strcat(folder,'/StaticParams.mat'),'-struct','StaticParams');

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
solveModel_NoAggShock(StaticParams,folder);

%--------------------------------------------------------------------------
% Generate orthogonal polynomials for polynomial chaos approximation
%--------------------------------------------------------------------------
% first set the order of truncation of the polynomial chaos approximation
order = 4;
setPoly(StaticParams,folder,order);

%--------------------------------------------------------------------------
% Model with aggregate risk
%--------------------------------------------------------------------------
solveModel_AggShock(StaticParams,case_nr,folder,order);

%--------------------------------------------------------------------------
% Compute Results
%--------------------------------------------------------------------------
results_statDistr(StaticParams,folder);
results_Errors(StaticParams,folder);