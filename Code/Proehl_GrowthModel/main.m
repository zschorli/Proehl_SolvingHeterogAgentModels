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
% 1 benchmark model (see online appendix on model parameters)
% 2 unemployment benefit of 40%
% 3 unemployment benefit of 65%

for case_nr = 1:3
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
order = 3;
setPoly(StaticParams,folder,order);

%--------------------------------------------------------------------------
% Model with aggregate risk
%--------------------------------------------------------------------------
solveModel_AggShock(StaticParams,case_nr,folder,order);
end

%--------------------------------------------------------------------------
% Run the Krusell-Smith algo for comparison
%--------------------------------------------------------------------------
for nr_moments=1:4
    file = strcat('C:\Users\Elisabeth Proehl\Documents\GitHub\Proehl_SolvingHeterogAgentModels\Code\KrusellSmithByMaliarMaliarValli_',num2str(nr_moments),'mom\MAIN.m');
    run(file);
end

%--------------------------------------------------------------------------
% Compute Results
%--------------------------------------------------------------------------
for case_nr = 1:3
folder = strcat('res_case',num2str(case_nr));        
StaticParams = load(strcat('C:\Users\Elisabeth Proehl\Documents\GitHub\Proehl_SolvingHeterogAgentModels\Code\Proehl_GrowthModel\res_case',num2str(case_nr),'\StaticParams.mat'));

results_statDistr(case_nr,StaticParams,folder);
results_Errors(case_nr,StaticParams,folder);
end
compareDistr_NoAggShock;