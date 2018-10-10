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
function [c_new,k_prime_new,y_new,y_new2,p_new,diffs,flag,sol_y,sol_y_2,sol_y_3,sol_nu,sol_nu_2,sol_nu_3,A2] ...
= OptPolUpdate_AggShock_PPA(pr,iter,agK_m,grid,...
  y,y2,sol_y,sol_y_1,sol_y_2,sol_y_3,sol_nu,sol_nu_2,sol_nu_3,A,StaticParams,Sol,Poly,c_pr_der)

diffs = ones(1,2);
b = 0.09;
k_prime = reshape(sol_y',[size(sol_y',1)/4,4])';

mult = 1.01;
m = mult^iter;
e = 1/(iter^1.51);
tol = (e^2)/(2*m); 
alpha = sqrt((b*A*m)^2+4*b*A*m)/2-b*A*m/2;

if iter<=15
    sol_k_large = interp1(StaticParams.kGrid_pol,k_prime',StaticParams.kGrid)';   
    [~,weights_new] = calcAggregates(sol_k_large,grid,StaticParams,Poly);
    weights_new(:,:,1) = 0;
    c_prime = @(z_idx,z,k,w)...
               interp1(StaticParams.kGrid_pol,...
               getCApprox(Sol.distrGrid,squeeze(sol_y_1(:,2*z(1)+z(2)+1,:)),...
               Sol.idx,squeeze(w((z_idx>2)+1,z(1)+1,:))'),...
               max(StaticParams.k_min,k),'linear','extrap');
    c_pr = zeros([4,size(sol_y_1,2),size(sol_y_1,3)]);
    for z_idx=1:4
        c_pr(z_idx,:,:) = [c_prime(z_idx,[0,0],k_prime(z_idx,:),weights_new);...
                           c_prime(z_idx,[0,1],k_prime(z_idx,:),weights_new);...
                           c_prime(z_idx,[1,0],k_prime(z_idx,:),weights_new);...
                           c_prime(z_idx,[1,1],k_prime(z_idx,:),weights_new)] ;
    end   
    sol_c_large = @(z_idx)...
                  interp1(StaticParams.kGrid_pol,...
                  (StaticParams.beta*(StaticParams.P(z_idx,:)*...
                  (squeeze(c_pr(z_idx,:,:)).^(-StaticParams.gamma)))).^(-1/StaticParams.gamma)...
                  ,StaticParams.kGrid);   
    C = calcAggregates([sol_c_large(1);sol_c_large(2);sol_c_large(3);sol_c_large(4)],grid,StaticParams,Poly);
    E = [squeeze(Sol.wage(agK_m,1:2,1))*StaticParams.p(1:2)*2;...
         squeeze(Sol.wage(agK_m,3:4,1))*StaticParams.p(3:4)*2];
    p_start = (E'./C').^StaticParams.gamma;
else
    p_start = pr;
end

if iter<=5
    p_new = p_start;
    [~,flag,sol] = solveEquCond(tol,k_prime,sol_y_1,p_new,y,y2,...
                      agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der,sol_y,iter);
else
    options = optimoptions('fsolve','Display','off','FunctionTolerance',min(StaticParams.criter_pdf,tol^2),'TolX',eps,'MaxIter',400);
    p_new = fsolve(@(p) solveEquCond(tol,k_prime,sol_y_1,p,y,y2,agK_m,...
                        grid,m,StaticParams,Sol,Poly,c_pr_der,sol_y,iter),...
                   p_start,options);
    [~,flag,sol] = solveEquCond(tol,k_prime,sol_y_1,p_new,y,y2,...
                          agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der,sol_y,iter);
    if flag <=0 && iter<=20
        p_new = p_start;
        [~,flag,sol] = solveEquCond(tol,k_prime,sol_y_1,p_new,y,y2,...
                          agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der,sol_y,iter);
        flag = 100+flag;
    end
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

function [A,flag,sol] = solveEquCond(tol,k_prime,sol_y_1,p,y,y2,agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der,sol_y)
   options = optimset('Display','off','Jacobian','on','TolFun',(tol^2),'TolX',eps,'MaxIter',400);
   [sol,~,flag] = fsolve(@(sol) FOC_AggShock(sol,k_prime,...
                          sol_y_1,p,y,y2,agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der),...
                          sol_y,options);
   
   k_prime_new = reshape(sol,size(k_prime'))';
   k_prime_large = interp1(StaticParams.kGrid_pol,k_prime_new',StaticParams.kGrid)';   
   A = calcAggregates(k_prime_large,grid,StaticParams,Poly);
end