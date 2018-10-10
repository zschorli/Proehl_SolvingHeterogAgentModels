% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This file implements the policy function iteration for the model without
% aggregate shocks.
%__________________________________________________________________________
function [c,k_prime,diffs] ...
         = OptPol_NoAggShock_PFI(bool_gb,c,k_prime,StaticParams,preSol)

grid = linspace(eps,100,10001);    
u_der_grid = 1./(grid.^StaticParams.gamma);         
u_der = @(x) interp1(grid,u_der_grid,max(eps,x),'linear','extrap');
u_der_inv = @(x) interp1(u_der_grid,grid,max(eps,x),'linear','extrap');
diffs = 1e8.*ones(2,1);
iter=0;
c_new = zeros(size(c));
if bool_gb
   P = StaticParams.P_g;
else
   P = StaticParams.P_b;
end   
	
while diffs(max(iter,1)) > StaticParams.criter_k*1e-3
   c_new(1,:,:) = u_der_inv(StaticParams.beta...
                ./repmat(permute(preSol.KGrid,[3,1,2]),[1,length(StaticParams.kGrid_pol),1])...
				.*sum(repmat(P(1,:)',[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)])...
                .*u_der(fct_k(StaticParams.kGrid_pol,c,k_prime([1 1],:,:))),1));
   c_new(2,:,:) = u_der_inv(StaticParams.beta...
                /repmat(permute(preSol.KGrid,[3,1,2]),[1,length(StaticParams.kGrid_pol),1])...
				.*sum(repmat(P(2,:)',[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)])...
				.*u_der(fct_k(StaticParams.kGrid_pol,c,k_prime([2 2],:,:))),1));
   c_new = max(0,c_new);
   k_prime_new = max(StaticParams.k_min,(preSol.wealth(2*bool_gb+(1:2),:,:)-c_new)...
                 ./repmat(permute(preSol.KGrid,[3,1,2]),[2,length(StaticParams.kGrid_pol),1]));
   c_new = preSol.wealth(2*bool_gb+(1:2),:,:)-k_prime_new.*repmat(permute(preSol.KGrid,[3,1,2]),[1,length(StaticParams.kGrid_pol),1]);
             
   iter = iter+1;
   diffs(iter) = max([max(max(squeeze(abs(k_prime_new(1,:,:)-k_prime(1,:,:))),[],2)),...
			          max(max(squeeze(abs(k_prime_new(2,:,:)-k_prime(2,:,:))),[],2))]);
   c = 0.5*c +0.5*c_new;
   k_prime = 0.5*k_prime +0.5*k_prime_new;
end
end

