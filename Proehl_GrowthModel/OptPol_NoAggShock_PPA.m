% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This file implements the proximal point algorithm for the model without
% aggregate shocks.
%__________________________________________________________________________
function [c,k_prime,y,diffs,flag] ...
         = OptPol_NoAggShock_PPA(i,bool_gb,c,k_prime,y,StaticParams,preSol,iter_ct)

sol_x = reshape(k_prime',1,numel(k_prime))';
sol_nu = sol_x;
sol_x_2 = reshape(y',1,numel(y))';
sol_nu_2 = sol_x_2;
A = 10;
b = 0.09;

iter = 1;
diffs = [1e8*ones(1,2);zeros(iter_ct-1,2)];
flag = zeros(iter_ct,1);
            
while max(diffs(max(1,iter-1),:)) > StaticParams.criter_k && iter-1<iter_ct     
     mult = 1.07;
     m = mult^iter;
     e = 1/(iter^1.51);
     tol = (e^2)/(2*m); 
     alpha = sqrt((b*A*m)^2+4*b*A*m)/2-b*A*m/2;
     sol_x = (1-alpha)*max(StaticParams.k_min,reshape(k_prime',1,numel(k_prime))')...
             +alpha*max(StaticParams.k_min,sol_nu);
     sol_x_1 = squeeze(preSol.wealth(bool_gb*2+(1:2),:,i)) - reshape(sol_x,size(k_prime'))';
     sol_x_2 = (1-alpha)*reshape(y',1,numel(y))'+alpha*sol_nu_2;
     options = optimset('Display','off','Jacobian','on','TolFun',tol^2,'TolX',eps,'MaxIter',500);
     [sol,~,flag(iter,:)] = fsolve(@(sol) FOC_NoAggShock(sol,sol_x_1,reshape(sol_x,size(k_prime'))',reshape(sol_x_2,size(y'))',i,bool_gb,m,StaticParams,preSol),...
                    sol_x,options);   
     k_prime_new = reshape(sol,size(k_prime'))';
     c_new = squeeze(preSol.wealth(bool_gb*2+(1:2),:,i)) - k_prime_new;
     y_new = max(0,y-m.*k_prime_new);
     
     sol_nu = sol_nu - (alpha/((1-alpha)*A*m))*(sol_x-sol);
     sol_nu_2 = sol_nu_2 - (alpha/((1-alpha)*A*m))*(sol_x_2-reshape(y_new',1,numel(y_new))');
     A = (1-alpha)*A;
     
     diffs(iter,1) = max(max(abs(k_prime_new-k_prime),[],2));
     diffs(iter,2) = max(max(abs(y_new-y),[],2));

     k_prime = k_prime_new;
     c = c_new;
     y = y_new;
     
     iter = iter+1;
end
end

