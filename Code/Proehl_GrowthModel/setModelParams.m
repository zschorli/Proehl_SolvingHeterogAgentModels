% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file sets the model parameters.
%__________________________________________________________________________

function StaticParams = setModelParams(StaticParams,case_nr,mu)
%--------------------------------------------------------------------------
% Model parameters given by den Haan, Judd and Juillard (2010)
%--------------------------------------------------------------------------
% unemployment benefit rate
if nargin>2
    StaticParams.mu = mu;
else
switch case_nr
    case 2
        StaticParams.mu = 0.4;
    case 3
        StaticParams.mu = 0.65;
    otherwise 
        StaticParams.mu = 0.15;
end
end

% risk aversion parameter
% switch case_nr
%     case 7
%         StaticParams.gamma=2;    
%     case 8
%         StaticParams.gamma=4;    
%     otherwise
        StaticParams.gamma=1;       
% end


StaticParams.beta=0.99; % time discount factor in the utility
StaticParams.alpha=0.36; % share of capital in the production function
StaticParams.delta=0.025; % depreciation rate
StaticParams.l_bar=1/0.9;  % time factor for labor
StaticParams.delta_a = 0.01; % change of aggregate productivity
StaticParams.kss_b = ((1/StaticParams.beta-(1-StaticParams.delta))/(StaticParams.alpha.*(1-StaticParams.delta_a)))^(1/(StaticParams.alpha-1));
StaticParams.kss_g = ((1/StaticParams.beta-(1-StaticParams.delta))/(StaticParams.alpha.*(1+StaticParams.delta_a)))^(1/(StaticParams.alpha-1));

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
% Individual capital grid
%--------------------------------------------------------------------------                
StaticParams.k_min = 0;          % minimum grid-value of individual capital

switch case_nr
    case 1
        StaticParams.k_max = 500;        % maximum grid-value of individual capital
        % grid for the policy      
        xgrid=linspace(0.5,15,45);      
        ygrid=(1-log(xgrid/xgrid(1))/log(xgrid(end)/xgrid(1))).^1.5; 
        ygrid=ygrid/max(ygrid); 
    case 2
        StaticParams.k_max = 800;        % maximum grid-value of individual capital
        % grid for the policy      
        xgrid=linspace(0.5,15,47);      
        ygrid=(1-log(xgrid/xgrid(1))/log(xgrid(end)/xgrid(1))).^1.5; 
        ygrid=ygrid/max(ygrid); 
    case 3
        StaticParams.k_max = 1000; 
        % grid for the policy      
        xgrid=linspace(0.5,15,49);      
        ygrid=(1-log(xgrid/xgrid(1))/log(xgrid(end)/xgrid(1))).^1.5; 
        ygrid=ygrid/max(ygrid); 
end
   
bound = StaticParams.k_min+2;
StaticParams.kGrid_pol = unique(sort([linspace(StaticParams.k_min,bound,36),...
                                 bound+(StaticParams.k_max-bound)*ygrid,...
                                 StaticParams.kss_b,StaticParams.kss_g]));   
% grid for the distribution
StaticParams.kGrid = linspace(StaticParams.k_min,StaticParams.k_max,4501);
        
clear xgrid ygrid;

%--------------------------------------------------------------------------
% Aggregate capital grid
%--------------------------------------------------------------------------
StaticParams.KGrid = unique(sort([linspace(33,43,40),StaticParams.kss_b,StaticParams.kss_g])); % aggregate capital grid
                  
%--------------------------------------------------------------------------
% Parameters for the fixed point iteration
%--------------------------------------------------------------------------
StaticParams.criter_k = 5e-5;  % convergence criterion for the capital policy
StaticParams.criter_pdf = 5e-6;  % convergence criterion for the distributions

%--------------------------------------------------------------------------
% Productivities and Cash on hand
%--------------------------------------------------------------------------
% dimensions: z^ag x z^id (2x2=4), k grid, K grid
% interest rate 
StaticParams.rate = @(z_ag,K)...
              StaticParams.alpha.*...
              repmat((1+(-1+2*z_ag).*StaticParams.delta_a),[2,1,numel(K)])...
              .*(repmat(permute(K,[1,3,2]),[2,1,1])...
              ./(StaticParams.l_bar.*repmat(StaticParams.er(1+z_ag),[2,1,numel(K)])))...
              .^(StaticParams.alpha-1);   
% wage 
StaticParams.wage = @(z_ag,K)...
              (1-StaticParams.alpha).*...
              repmat((1+(-1+2*z_ag).*StaticParams.delta_a),[2,1,numel(K)])...
              .*(repmat(permute(K,[1,3,2]),[2,1,1])...
              ./(StaticParams.l_bar.*repmat(StaticParams.er(1+z_ag),[2,1,numel(K)])))...
              .^(StaticParams.alpha);   
% income
StaticParams.wealth = @(z_ag,K,k)...
              ((StaticParams.l_bar-...
                repmat(StaticParams.mu.*StaticParams.ur(1+z_ag)/StaticParams.er(1+z_ag),[2,numel(k),numel(K)]))...
                .*repmat([0;1],[1,numel(k),numel(K)])...
                +StaticParams.mu.*(1-repmat([0;1],[1,numel(k),numel(K)])))...
                .*repmat(StaticParams.wage(z_ag,K),[1,numel(k),1])...
                +(1-StaticParams.delta+repmat(StaticParams.rate(z_ag,K),[1,numel(k),1]))...
                .*repmat(k,[2,1,numel(K)]);
end