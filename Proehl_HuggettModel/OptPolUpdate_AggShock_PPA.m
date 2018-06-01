% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This file implements the update of the proximal point algorithm for the 
% model with aggregate shocks.
%__________________________________________________________________________
function [c_new,k_prime_new,y_new,y_new2,p_new,diffs,flag,sol_y,sol_y_2,sol_y_3,sol_nu,sol_nu_2,sol_nu_3,A2] ...
= OptPolUpdate_AggShock_PPA(pr,iter,agK_m,grid,...
  y,y2,sol_y,sol_y_1,sol_y_2,sol_y_3,sol_nu,sol_nu_2,sol_nu_3,A,StaticParams,Sol,Poly,c_pr_der)

diffs = ones(1,2);
b = 0.1;
k_prime = reshape(sol_y',[size(sol_y',1)/4,4])';

mult = 1.01;
m = mult^iter;
e = 1/(iter^2);
tol = (e^2)/(2*m); 
alpha = sqrt((b*A*m)^2+4*b*A*m)/2-b*A*m/2;

% outer loop for finding the price
if (Sol.diff(max(1,iter-1),3)>1e-10 || iter<=10) 
    options = optimset('Display','off','Jacobian','off','TolFun',min(1e-6,tol^2),'TolX',eps,'MaxIter',400);
    [p_new,~,flag_p] = fsolve(@(p) EC_AggShock(p,k_prime,sol_y,sol_y_1,y,y2,agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der,tol),...
                   pr,options);
else
    p_new = pr;
    flag_p = 1; 
end
options = optimset('Display','off','Jacobian','on','TolFun',(tol^2),'TolX',eps,'MaxIter',400);
[sol,~,flag_k] = fsolve(@(sol) FOC_AggShock(p_new,sol,k_prime,...
                      sol_y_1,y,y2,agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der),...
                      sol_y,options);
if min(flag_p,flag_k)<=0
   flag = min(flag_p,flag_k);                  
else
   flag = max(flag_p,flag_k);    
end

% dampening because initital iterations are too far from solution
if iter<=20
    p_new = 0.6*p_new+0.4*pr;
end
           
k_prime_new = reshape(sol,size(k_prime'))';
c_new = max(eps,squeeze(Sol.wealth(agK_m,:,:))-repmat(p_new([1 1 2 2])',[1,size(k_prime_new,2)]).*k_prime_new);
y_new = max(0,y-m.*(k_prime_new-StaticParams.k_min));
y_new2 = max(0,y2-m.*(squeeze(Sol.wealth(agK_m,:,:))-k_prime_new.*repmat(p_new([1 1 2 2])',[1,size(k_prime,2)])));

sol_nu = sol_nu - (alpha/((1-alpha)*A*m))*(sol_y-sol); 
sol_nu_2 = sol_nu_2 - (alpha/((1-alpha)*A*m))*(sol_y_2-reshape(y_new',1,numel(y_new)));      
sol_nu_3 = sol_nu_3 - (alpha/((1-alpha)*A*m))*(sol_y_3-reshape(y_new2',1,numel(y_new2)));      
A2 = (1-alpha)*A;
m = mult^(iter+1);
alpha = sqrt((b*A2*m)^2+4*b*A2*m)/2-b*A2*m/2;
sol_y = (1-alpha)*max(StaticParams.k_min,sol)...
             +alpha*max(StaticParams.k_min,sol_nu);
sol_y_2 = (1-alpha)*reshape(y_new',1,numel(y_new))...
               +alpha*sol_nu_2;
sol_y_3 = (1-alpha)*reshape(y_new2',1,numel(y_new2))...
               +alpha*sol_nu_3;  

diffs(1) = max(max(max(abs(k_prime_new-k_prime))));
diffs(2) = max(max(max(abs(y_new-y))));
diffs(3) = max(max(max(abs(p_new-pr))));
end

function S = EC_AggShock(p,k_prime,sol_y,sol_y_1,y,y2,agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der,tol)
% inner loop for finding the policy
options = optimset('Display','off','Jacobian','on','TolFun',(tol^2),'TolX',eps,'MaxIter',400);
sol = fsolve(@(sol) FOC_AggShock(p,sol,k_prime,...
             sol_y_1,y,y2,agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der),...
             sol_y,options);  
         
k_prime_new = reshape(sol,size(k_prime'))';
k_prime_large = interp1(StaticParams.kGrid_pol,k_prime_new',StaticParams.kGrid)';
S = calcAggregates(k_prime_large,grid,StaticParams,Poly);
end