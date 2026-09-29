% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This file sets the model parameters.
%__________________________________________________________________________

function StaticParams = setModelParams(StaticParams)
%--------------------------------------------------------------------------
% Model parameters given by den Krusell, Mukoyama, and Smith (2011)
%--------------------------------------------------------------------------
% stationary distribution of exogenous shocks       
StaticParams.z=[0 0;0 1;1 0;1 1]; % aggregate and idiosyncratic shocks  
StaticParams.ur = [0.5;0.5]; 
StaticParams.er = 1-StaticParams.ur;     
StaticParams.p = [[StaticParams.ur(1);StaticParams.er(1)]*0.5;...
                  [StaticParams.ur(2);StaticParams.er(2)]*0.5];  

% transition probability matrix (states: bad-unempl, bad-empl, good-unempl,
% good-empl)          
StaticParams.P = [0.942 0.058 0.984 0.016;...
                  0.058 0.942 0.016 0.984;...
                  0.942 0.058 0.984 0.016;...
                  0.058 0.942 0.016 0.984];
StaticParams.P = StaticParams.P.*[0.43*ones(2,2),0.57*ones(2,2);...
                                  0.57*ones(2,2),0.43*ones(2,2)];
StaticParams.P_g = StaticParams.P(3:4,3:4)./repmat(sum(StaticParams.P(3:4,3:4),2),[1,2]);
StaticParams.P_b = StaticParams.P(1:2,1:2)./repmat(sum(StaticParams.P(1:2,1:2),2),[1,2]);

%--------------------------------------------------------------------------
% Model parameters given by Kehoe, Midrigan, and Pastorino (2016)
%--------------------------------------------------------------------------
% risk aversion parameter
StaticParams.gamma = 2;    
% elasticity of substitution
StaticParams.eps = 1./StaticParams.gamma;
% collateral constraint
StaticParams.theta = 0.8;

%--------------------------------------------------------------------------
% Model parameters calibrated to match 2% interest rate and a 
% housing-income-ratio of 2.2 (when taking the average of the two models 
% without aggregat risk, i.e. fixing good/bad state)
%--------------------------------------------------------------------------
% time discount factor in the utility
StaticParams.beta= 0.964;
% preference parameter for housing
StaticParams.eta = 0.248;

%--------------------------------------------------------------------------
% Individual bond holdings grid
%--------------------------------------------------------------------------      
% minimum grid-value of individual asset holdings
StaticParams.a_min = 0; 

% maximum grid-value of individual asset holdings
StaticParams.a_max = 60;       

% grid for the policy    
xgrid=linspace(0.5,15,50);
ygrid=(1-log(xgrid/xgrid(1))/log(xgrid(end)/xgrid(1))).^1.5; 
ygrid=ygrid/max(ygrid);    
bound = 1;%StaticParams.a_min+0.5*abs(StaticParams.a_min);
StaticParams.aGrid_pol = unique([linspace(StaticParams.a_min,bound,36),...
                                 bound+(StaticParams.a_max-bound)*ygrid]);
StaticParams.aGrid_pol(end) = StaticParams.a_max;                          

% grid for the distribution
StaticParams.aGrid = unique(linspace(StaticParams.a_min,StaticParams.a_max,1001));

StaticParams.hGrid = linspace(0,2.5,251);
StaticParams.cGrid = linspace(0.4,2.4,251);
StaticParams.bGrid = linspace(-4,14,251);

clear xgrid ygrid;
               
%--------------------------------------------------------------------------
% Aggregate grids
%--------------------------------------------------------------------------
StaticParams.NA = 101;
StaticParams.Np = 21;
StaticParams.RGrid = linspace(1+eps,1/0.963-eps,StaticParams.Np); %R=1+r
StaticParams.qGrid = linspace(1.8,3,StaticParams.NA); 
StaticParams.R_large = repmat(permute(StaticParams.RGrid,[3,1,2]),[2,length(StaticParams.aGrid_pol),1,StaticParams.NA]); 
StaticParams.q_large = repmat(permute(StaticParams.qGrid,[3,4,1,2]),[2,length(StaticParams.aGrid_pol),StaticParams.Np,1]); 
StaticParams.uq_large = (StaticParams.R_large-1).*StaticParams.q_large; 

%--------------------------------------------------------------------------
% Fixed quantities
%--------------------------------------------------------------------------
% endowment
StaticParams.wage = repmat([0.982*[2/3.06;4.12/3.06];...
                      1.054*[2/3.06;4.12/3.06]],[1,length(StaticParams.aGrid_pol)]);
% cash on hand
StaticParams.wealth = @(z_ag) repmat(StaticParams.wage(2*z_ag+(1:2),:),[1,1,StaticParams.Np,StaticParams.NA])...
                       +StaticParams.R_large.*repmat(StaticParams.aGrid_pol,[2,1,StaticParams.Np,StaticParams.NA]);

%--------------------------------------------------------------------------
% Model parameters for the fixed point iteration
%--------------------------------------------------------------------------
StaticParams.criter_a = 5e-4;  % convergence criterion for the capital policy
StaticParams.criter_pdf = 5e-7;  % convergence criterion for the distributions
end