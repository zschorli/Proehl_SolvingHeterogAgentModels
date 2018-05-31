% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This is the main file
%__________________________________________________________________________

%--------------------------------------------------------------------------
% Choose model configuration
%--------------------------------------------------------------------------
% Choose a number between 1 and 20:
% 1 benchmark model (see online appendix on model parameters)
% 2 unemployment benefit of 65%
% 3 unemployment benefit of 5%
% 4 unemployment benefit of 40%
% 5 subjective discount factor of 0.9
% 6 subjective discount factor of 0.999
% 7 risk aversion of 2
% 8 risk aversion of 4
% 9 production function elasticity of 0.2
% 10 production function elasticity of 0.45
% 11 depreciation rate of 1%
% 12 depreciation rate of 7.5%
% 13 unemployment benefit of 65%, subjective discount factor of 0.9
% 14 unemployment benefit of 65%, subjective discount factor of 0.999
% 15 unemployment benefit of 65%, risk aversion of 2
% 16 unemployment benefit of 65%, risk aversion of 4
% 17 unemployment benefit of 65%, production function elasticity of 0.2
% 18 unemployment benefit of 65%, production function elasticity of 0.45
% 19 unemployment benefit of 65%, depreciation rate of 1%
% 20 unemployment benefit of 65%, depreciation rate of 7.5%

case_nr = 1;

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
solveModel_NoAggShock(StaticParams,folder,case_nr);

%--------------------------------------------------------------------------
% Generate orthogonal polynomials for polynomial chaos approximation
%--------------------------------------------------------------------------
% first set the order of truncation of the polynomial chaos approximation
switch case_nr
    case {2,4}        
        order = 4;
    otherwise 
        order = 3;
end
setPoly(StaticParams,folder,order);

%--------------------------------------------------------------------------
% Model with aggregate risk
%--------------------------------------------------------------------------
solveModel_AggShock(StaticParams,case_nr,folder,order);

%--------------------------------------------------------------------------
% Compute Results
%--------------------------------------------------------------------------
results_statDistr(case_nr,StaticParams,folder);
if case_nr<=2
    results_Errors(case_nr,StaticParams,folder);
end
