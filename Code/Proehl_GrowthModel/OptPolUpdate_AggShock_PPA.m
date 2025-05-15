% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file implements the update of the proximal point algorithm for the 
% model with aggregate shocks.
%__________________________________________________________________________
function [c_new,k_prime_new,y_new,diffs,flag] ...
= OptPolUpdate_AggShock_PPA(iter,agK_m,grid,k_prime,y,...
  StaticParams,Sol,Poly,m,agK_prime_new_dk,c_pr_der)

diffs = ones(1,2);
% if M<=1
%     m=0.1;
% else
%     m = 0.2;%min(0.4,0.1+((iter-1)*5e-3));
% end

update = min(1,0.05*(1.5^(iter-1)));
      
options = optimoptions(@fsolve,'Display','off','SpecifyObjectiveGradient',true,...
          'MaxIterations',100,'FunctionTolerance',1e-6,'StepTolerance',1e-10);
sol_k = reshape(k_prime',[],1);
[sol,~,flag] = fsolve(@(sol) FOC_AggShock(sol,k_prime,y,agK_m,grid,m,...
                             StaticParams,Sol,Poly,agK_prime_new_dk,c_pr_der),...
                             sol_k,options);	
k_prime_aux = reshape(sol,size(k_prime'))';
y_aux = max(0,y-m.*k_prime_aux);
diffs(1) = max(max(abs(k_prime_aux-k_prime),[],2));
diffs(2) = max(max(abs(y_aux-y),[],2));

k_prime_new = update.*k_prime_aux+(1-update).*k_prime;
c_new = squeeze(Sol.wealth(agK_m,:,:))-k_prime_new;
y_new = max(0,y-m.*k_prime_new);
end