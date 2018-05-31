% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This file sets the model parameters.
%__________________________________________________________________________

function StaticParams = setModelParams(StaticParams,case_nr)
%--------------------------------------------------------------------------
% Model parameters given by den Haan, Judd and Juillard (2010)
%--------------------------------------------------------------------------
% unemployment benefit rate
switch case_nr
    case {2,13,14,15,16,17,18,19,20}
        StaticParams.mu = 0.65;
    case 3
        StaticParams.mu = 0.05;
    case 4
        StaticParams.mu = 0.4;
    otherwise 
        StaticParams.mu = 0.15;
end

% time discount factor in the utility
switch case_nr
    case {5,13}
        StaticParams.beta=0.9; 
    case {6,14}
        StaticParams.beta=0.999; 
    otherwise
        StaticParams.beta=0.99;    
end

% risk aversion parameter
switch case_nr
    case {7,15}
        StaticParams.gamma=2;    
    case {8,16}
        StaticParams.gamma=4;    
    otherwise
        StaticParams.gamma=1;       
end

% share of capital in the production function
switch case_nr
    case {9,17}
        StaticParams.alpha=0.2;
    case {10,18}
        StaticParams.alpha=0.45;
    otherwise
        StaticParams.alpha=0.36;    
end

% depreciation rate
switch case_nr
    case {11,19}
        StaticParams.delta=0.01;
    case {12,20}
        StaticParams.delta=0.075;
    otherwise
        StaticParams.delta=0.025;
end
        
StaticParams.l_bar=1/0.9;  % time factor for labor
StaticParams.delta_a = 0.01; % change of aggregate productivity

% stationary distribution of exogenous shocks     
StaticParams.z=[0 0;0 1;1 0;1 1]; % aggregate and idiosyncratic shocks    
StaticParams.ur = [0.1;0.04]; 
StaticParams.er = 1-StaticParams.ur;     
StaticParams.p = [[StaticParams.ur(1);StaticParams.er(1)]*0.5;...
                  [StaticParams.ur(2);StaticParams.er(2)]*0.5];  

% transition probability matrix (states: bad-unempl, bad-empl, good-unempl,
% good-empl)          
StaticParams.P = [0.525 0.35 0.03125 0.09375  
                  0.038889 0.836111 0.002083 0.122917
                  0.09375 0.03125 0.291667 0.583333
                  0.009115 0.115885 0.024306 0.850694];
StaticParams.P_g = StaticParams.P(3:4,3:4)./repmat(sum(StaticParams.P(3:4,3:4),2),[1,2]);
StaticParams.P_b = StaticParams.P(1:2,1:2)./repmat(sum(StaticParams.P(1:2,1:2),2),[1,2]);
              
%--------------------------------------------------------------------------
% Steady state capital
%--------------------------------------------------------------------------
% steady-state capital in a deterministic model according to Maliar, Maliar
% and Valli (2010)
StaticParams.kss=((1/StaticParams.beta-(1-StaticParams.delta))/StaticParams.alpha)^...
                (1/(StaticParams.alpha-1)); 

%--------------------------------------------------------------------------
% Individual capital grid
%--------------------------------------------------------------------------                
StaticParams.k_min = 0;          % minimum grid-value of individual capital
StaticParams.k_max = 450;        % maximum grid-value of individual capital

% grid for the policy      
xgrid=linspace(0.5,15,45);      
ygrid=(1-log(xgrid/xgrid(1))/log(xgrid(end)/xgrid(1))).^1.5; 
ygrid=ygrid/max(ygrid);    
bound = 2;
StaticParams.kGrid_pol = unique([linspace(StaticParams.k_min,bound,36),...
                                 bound+(StaticParams.k_max-bound)*ygrid]);   
% grid for the distribution
StaticParams.kGrid = unique(linspace(StaticParams.k_min,StaticParams.k_max,...
                            (StaticParams.k_max-StaticParams.k_min)*10+1));
        
clear xgrid ygrid;
                  
%--------------------------------------------------------------------------
% Parameters for the fixed point iteration
%--------------------------------------------------------------------------
StaticParams.methodEqVFI = 1; 
StaticParams.criter_k = 5e-5;  % convergence criterion for the capital policy
StaticParams.criter_pdf = 5e-7;  % convergence criterion for the distributions
end