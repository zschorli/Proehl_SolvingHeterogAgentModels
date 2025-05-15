% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file implements the proximal point algorithm for the model without
% aggregate shocks.
%__________________________________________________________________________
function [c,k_prime,y,diffs,flag] ...
         = OptPol_NoAggShock_PPA(i,bool_gb,k_prime,y,StaticParams,iter_ct)

iter = 1;
diffs = [1e8*ones(1,2);zeros(iter_ct-1,2)];
flag = zeros(iter_ct,1);
m = 0.5;
update = 0.05;
options = optimoptions(@fsolve,'Display','off','SpecifyObjectiveGradient',true,...
          'MaxIterations',100,'FunctionTolerance',1e-6,'StepTolerance',1e-10);

while max(diffs(max(1,iter-1),:)) > StaticParams.criter_k && iter-1<iter_ct 
     [sol,~,flag(iter,:)] = fsolve(@(sol) FOC_NoAggShock(sol,k_prime,y,i,bool_gb,m,StaticParams),...
                    reshape(k_prime,[],1),options);  
     % display(flag(iter,:));
     k_prime_new = reshape(sol,size(k_prime));
     y_new = max(0,y-m.*k_prime_new);
          
     diffs(iter,1) = max(max(abs(k_prime_new-k_prime),[],2));
     diffs(iter,2) = max(max(abs(y_new-y),[],2));

     k_prime = update*k_prime_new+(1-update)*k_prime;
     y = max(0,y-m.*k_prime);
     
     update = min(1,update*1.5);
     m = min(3,m+5e-4);
     iter = iter+1;
end
k_prime = max(0,k_prime);
c = StaticParams.wealth(bool_gb,StaticParams.KGrid(i),StaticParams.kGrid_pol) - k_prime;
end

