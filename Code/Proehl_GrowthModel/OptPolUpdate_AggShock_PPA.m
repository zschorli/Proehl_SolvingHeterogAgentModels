% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This file implements the update of the proximal point algorithm for the 
% model with aggregate shocks.
%__________________________________________________________________________
function [c_new,k_prime_new,y_new,diffs,flag,sol_y,sol_y_2,sol_nu,sol_nu_2,A] ...
= OptPolUpdate_AggShock_PPA(iter,agK_m,grid,...
  k_prime,y,sol_y,sol_y_1,sol_y_2,sol_nu,sol_nu_2,A,StaticParams,Sol,Poly,agK_prime_new_dk,c_pr_der)

diffs = ones(1,2);
b = 0.09;

mult = 1.07;
m = mult^iter;
e = 1/(iter^1.51);
tol = (e^2)/(2*m); 
alpha = sqrt((b*A*m)^2+4*b*A*m)/2-b*A*m/2;
      
options = optimset('Display','off','Jacobian','on','TolFun',tol^2,'TolX',eps,'MaxIter',400);
[sol,~,flag] = fsolve(@(sol) FOC_AggShock(sol,reshape(sol_y,size(k_prime'))',sol_y_1,reshape(sol_y_2,size(y'))',agK_m,grid,m,StaticParams,Sol,Poly,agK_prime_new_dk,c_pr_der),...
               sol_y,options);	
k_prime_new = reshape(sol,size(k_prime'))';
c_new = squeeze(Sol.wealth(agK_m,:,:))-k_prime_new;
y_new = max(0,y-m.*k_prime_new);

sol_nu = sol_nu - (alpha/((1-alpha)*A*m))*(sol_y-sol);
sol_nu_2 = sol_nu_2 - (alpha/((1-alpha)*A*m))*(sol_y_2-reshape(y_new',numel(y_new),1));      
A = (1-alpha)*A;
m = mult^(iter+1);
alpha = sqrt((b*A*m)^2+4*b*A*m)/2-b*A*m/2;
sol_y = (1-alpha)*max(StaticParams.k_min,reshape(k_prime_new',numel(k_prime_new),1))...
        +alpha*max(StaticParams.k_min,sol_nu);
sol_y_2 = (1-alpha)*reshape(y_new',numel(y_new),1)...
          +alpha*sol_nu_2;

diffs(1) = max(max(abs(k_prime_new-k_prime),[],2));
diffs(2) = max(max(abs(y_new-y),[],2));
end