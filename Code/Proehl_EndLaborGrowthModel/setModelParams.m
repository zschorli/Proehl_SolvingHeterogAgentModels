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
        case {1,2,14,19,21} 
            StaticParams.mu = 0.15;
        case {3,4}
            StaticParams.mu = 0.4;
        case 5
            StaticParams.mu = 0.135;
        case 6
            StaticParams.mu = 0.165;
        case {7,9,11,13,15,17,20}
            StaticParams.mu = 0;
        case {8,12,16,18}
            StaticParams.mu = 0.015;
        case 10
            StaticParams.mu = 0.03; 
    end
end
switch case_nr
    case {5,9}
        StaticParams.mu_extra = 0.03;
    case {6,10}
        StaticParams.mu_extra = -0.03;
    case {11,17}
        StaticParams.mu_extra = 0.015;
    case {12,18}
        StaticParams.mu_extra = -0.015;
    case {13,20}
        StaticParams.mu_extra = 0.15;
    case {14,21}
        StaticParams.mu_extra = -0.15;
    otherwise
        StaticParams.mu_extra = 0;
end

StaticParams.gamma=1;       
StaticParams.theta = 1/2.9; % Cobb-Douglas utility: c^theta*(1-l)^(1-theta)
StaticParams.beta=0.99; % time discount factor in the utility
StaticParams.alpha=0.36; % share of capital in the production function
StaticParams.delta=0.025; % depreciation rate
StaticParams.l_bar=1/0.9;  % time factor for labor
StaticParams.delta_a = 0.01; % change of aggregate productivity
switch case_nr
    case {2,4}
        StaticParams.rho = 0; %tax progressivity level (linear tax if rho=0)
    otherwise
        StaticParams.rho = 0.151; 
end

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
% Individual grids
%--------------------------------------------------------------------------      
if case_nr<=14
    StaticParams.k_min = 0;          % minimum grid-value of individual capital
else
    StaticParams.k_min = 20;        
end
if StaticParams.rho == 0
    StaticParams.k_max = 100;
else
    StaticParams.k_max = 170;        % maximum grid-value of individual capital
end
% grid for the policy      
xgrid=linspace(0.5,15,47);      
ygrid=(1-log(xgrid/xgrid(1))/log(xgrid(end)/xgrid(1))).^1.5; 
ygrid=ygrid/max(ygrid); 
bound = StaticParams.k_min+0.2;
StaticParams.kGrid_pol = unique(sort([linspace(StaticParams.k_min,bound,36),...
                                 bound+(StaticParams.k_max-bound)*ygrid])); 

StaticParams.lu = 0.2;
StaticParams.l_min = 0.25;

% grid for the distributions
StaticParams.kGrid = linspace(StaticParams.k_min,StaticParams.k_max,4501);
StaticParams.lGrid = linspace(0,1,501);
StaticParams.cGrid = linspace(0,3,1001);
clear xgrid ygrid;

%--------------------------------------------------------------------------
% Aggregate grids
%--------------------------------------------------------------------------
StaticParams.KGrid = linspace(9.5,14.5,51); % aggregate capital grid
StaticParams.LGrid = linspace(0.2,0.7,11); % aggregate labor grid
if StaticParams.rho == 0
    StaticParams.LtaxGrid = 1;
else
    StaticParams.LtaxGrid = linspace(0.75,0.9,6); % aggregate labor/ aggregate labor^(1-rho) grid
end
                  
%--------------------------------------------------------------------------
% Parameters for the fixed point iteration
%--------------------------------------------------------------------------
StaticParams.criter_k = 5e-5;  % convergence criterion for the capital policy
StaticParams.criter_pdf = 5e-6;  % convergence criterion for the distributions

%--------------------------------------------------------------------------
% Productivities and Cash on hand
%--------------------------------------------------------------------------
% dimensions: z^ag x z^id (2x2=4), k grid, K grid, L grid
% interest rate 
StaticParams.rate = @(z_ag,K,L)...
              StaticParams.alpha.*...
              repmat((1+(-1+2*z_ag).*StaticParams.delta_a),[2,1,numel(K),numel(L)])...
              .*(repmat(permute(K,[3,1,2]),[2,1,1,numel(L)])...
              ./(StaticParams.l_bar.*repmat(StaticParams.er(1+z_ag),[2,1,numel(K),numel(L)])...
                 .*repmat(permute(L,[3,4,1,2]),[2,1,numel(K),1])))...
              .^(StaticParams.alpha-1);   
% wage 
StaticParams.wage = @(z_ag,K,L)...
              (1-StaticParams.alpha).*...
              repmat((1+(-1+2*z_ag).*StaticParams.delta_a),[2,1,numel(K),numel(L)])...
              .*(repmat(permute(K,[3,1,2]),[2,1,1,numel(L)])...
              ./(StaticParams.l_bar.*repmat(StaticParams.er(1+z_ag),[2,1,numel(K),numel(L)])...
                 .*repmat(permute(L,[3,4,1,2]),[2,1,numel(K),1])))...
              .^(StaticParams.alpha);  

% tax rate (1-tau)
% LtaxRatio = L/(int l^(1-rho) d mu) < 1
StaticParams.taxrate = @(z_ag,K,L,LtaxRatio)...
              (1-repmat((StaticParams.mu+(z_ag==0).*StaticParams.mu_extra).*StaticParams.ur(1+z_ag).*StaticParams.lu,[2,1,numel(K),numel(L),numel(LtaxRatio)])...
                ./repmat(permute(StaticParams.er(1+z_ag).*StaticParams.l_bar.*L,[3,4,1,2]),[2,1,numel(K),1,numel(LtaxRatio)]))...
                .*(StaticParams.l_bar.*repmat(StaticParams.wage(z_ag,K,L),[1,1,1,1,numel(LtaxRatio)])).^StaticParams.rho...
                .*repmat(permute(LtaxRatio,[3,4,5,1,2]),[2,1,numel(K),numel(L),1]);

% income
StaticParams.wealth = @(z_ag,K,L,LtaxRatio,k,l)...% l has the same size as k
              repmat(StaticParams.taxrate(z_ag,K,L,LtaxRatio),[1,numel(k),1,1,1,1])...
                .*repmat([0;1],[1,numel(k),numel(K),numel(L),numel(LtaxRatio)])...
                .*(StaticParams.l_bar.*repmat(StaticParams.wage(z_ag,K,L),[1,numel(k),1,1,numel(LtaxRatio)]).*l).^(1-StaticParams.rho)...
              +(StaticParams.mu+(z_ag==0).*StaticParams.mu_extra).*(1-repmat([0;1],[1,numel(k),numel(K),numel(L),numel(LtaxRatio)]))...
                .*repmat(StaticParams.wage(z_ag,K,L),[1,numel(k),1,1,numel(LtaxRatio)]).*StaticParams.lu...
              +(1-StaticParams.delta+repmat(StaticParams.rate(z_ag,K,L),[1,numel(k),1,1,numel(LtaxRatio)]))...
                .*repmat(k,[2,1,numel(K),numel(L),numel(LtaxRatio)]);

% FOC l for employed labor choice: 
% c = (1-l)*wealth_derl*(theta/(1-theta)) = (1-l)*l^(-rho)*Q
% such that c/(1-l)=l^(-rho)*Q, i.e., Cobb-Douglas utility is
% u = c^theta*(1-l)^(1-theta) = (1-l)*l^(-rho*theta)*Q^theta
StaticParams.Q = @(z_ag,K,L,LtaxRatio)...
              StaticParams.taxrate(z_ag,K,L,LtaxRatio)...
                .*repmat([0;1],[1,1,numel(K),numel(L),numel(LtaxRatio)])...
                .*(StaticParams.l_bar.*repmat(StaticParams.wage(z_ag,K,L),[1,1,1,1,numel(LtaxRatio)])).^(1-StaticParams.rho)...
              .*(1-StaticParams.rho)...
              .*(StaticParams.theta./(1-StaticParams.theta)); 

% output 
StaticParams.output = @(z_ag,K,L)...
              (1+(-1+2*z_ag).*StaticParams.delta_a)...
              .*(K).^(StaticParams.alpha)...
              .*(StaticParams.l_bar.*StaticParams.er(1+z_ag)...
                 .*L).^(1-StaticParams.alpha);
end