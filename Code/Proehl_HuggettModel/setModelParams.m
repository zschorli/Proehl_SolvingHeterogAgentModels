% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This file sets the model parameters.
%__________________________________________________________________________

function StaticParams = setModelParams(StaticParams,case_nr)
%--------------------------------------------------------------------------
% Model parameters given by den Krusell, Mukoyama, and Smith (2011)
%--------------------------------------------------------------------------
% time discount factor in the utility
switch case_nr
    case {6,10}
        StaticParams.beta=0.9; 
    case {7,11}
        StaticParams.beta=0.999; 
    otherwise
        StaticParams.beta=0.99;
end
% risk aversion parameter
switch case_nr
    case {8,12}
        StaticParams.gamma=1;    
    case {9,13}
        StaticParams.gamma=3;   
    otherwise
        StaticParams.gamma=2;      
end

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
% Individual bond holdings grid
%--------------------------------------------------------------------------      
% minimum grid-value of individual asset holdings
switch case_nr
    case {2,10,11,12,13}
        StaticParams.k_min = -0.5; 
    case 3
        StaticParams.k_min = -1; 
    case 4
        StaticParams.k_min = -2; 
    case 5
        StaticParams.k_min = -2.5; 
    otherwise
        StaticParams.k_min = -1.5; 
end

% maximum grid-value of individual asset holdings
StaticParams.k_max = 3*abs(StaticParams.k_min);         

% grid for the policy    
xgrid=linspace(0.5,15,45);
ygrid=(1-log(xgrid/xgrid(1))/log(xgrid(end)/xgrid(1))).^1.5; 
ygrid=ygrid/max(ygrid);    
bound = StaticParams.k_min+0.5*abs(StaticParams.k_min);
StaticParams.kGrid_pol = unique([linspace(StaticParams.k_min,bound,36),...
                                 bound+(StaticParams.k_max-bound)*ygrid]);
StaticParams.kGrid_pol(end) = StaticParams.k_max;                          

% grid for the distribution
StaticParams.kGrid = unique(linspace(StaticParams.k_min,StaticParams.k_max,151));%(StaticParams.k_max-StaticParams.k_min)*50+1));

clear xgrid ygrid;
                  
%--------------------------------------------------------------------------
% Model parameters for the fixed point iteration
%--------------------------------------------------------------------------
StaticParams.methodEqPFI = 1; 
StaticParams.criter_k = 5e-5;  % convergence criterion for the capital policy
StaticParams.criter_pdf = 5e-7;  % convergence criterion for the distributions
StaticParams.update_k = 0.7;   % updating parameter for the capital policy

end