% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This is the main file
%__________________________________________________________________________

folder = strcat('res_Housing');        
if exist(folder,'dir') ~= 7
    mkdir(folder);
end

%--------------------------------------------------------------------------
% Parameter definition
%--------------------------------------------------------------------------
StaticParams = struct;
StaticParams = setModelParams(StaticParams);
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
solveModel_AggShock(StaticParams,folder,order);

%--------------------------------------------------------------------------
% Compute Results
%--------------------------------------------------------------------------
results_statDistr(StaticParams,folder,order);
results_Errors(StaticParams,folder,order);

% plots;