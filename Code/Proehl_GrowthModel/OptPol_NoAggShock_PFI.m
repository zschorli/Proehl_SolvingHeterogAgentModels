% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file implements the policy function iteration for the model without
% aggregate shocks.
%__________________________________________________________________________
function [c,k_prime,diffs] ...
         = OptPol_NoAggShock_PFI(bool_gb,c,k_prime,StaticParams)
 
diffs = 1e8.*ones(2,1);
iter=0;
c_new = zeros(size(c));
if bool_gb
   P = StaticParams.P_g;
else
   P = StaticParams.P_b;
end   
MargUt = @(c) c.^(-StaticParams.gamma);
MargUt_inv = @(c) c.^(-1./StaticParams.gamma);
wealth = StaticParams.wealth(bool_gb,StaticParams.KGrid,StaticParams.kGrid_pol);
rate = repmat(StaticParams.rate(bool_gb,StaticParams.KGrid),[1,numel(StaticParams.kGrid_pol),1]);
	
while diffs(max(iter,1)) > StaticParams.criter_k/5
   c_new(1,:,:) = max(0,min(wealth(1,:,:)-StaticParams.k_min,MargUt_inv(StaticParams.beta...
				.*sum(repmat(P(1,:)',[1,length(StaticParams.kGrid_pol),length(StaticParams.KGrid)])...
                      .*(1-StaticParams.delta+rate(1,:,:))...
				      .*MargUt(fct_k(StaticParams.kGrid_pol,c,k_prime([1 1],:,:))),1))));
   c_new(2,:,:) = max(0,min(wealth(2,:,:)-StaticParams.k_min,MargUt_inv(StaticParams.beta...
				.*sum(repmat(P(2,:)',[1,length(StaticParams.kGrid_pol),length(StaticParams.KGrid)])...
				      .*(1-StaticParams.delta+rate(2,:,:))...
                      .*MargUt(fct_k(StaticParams.kGrid_pol,c,k_prime([2 2],:,:))),1))));
   k_prime_new = max(StaticParams.k_min,wealth-c_new);

   iter = iter+1;
   diffs(iter) = max([max(max(squeeze(abs(k_prime_new(1,:,:)-k_prime(1,:,:))),[],2)),...
			          max(max(squeeze(abs(k_prime_new(2,:,:)-k_prime(2,:,:))),[],2))]);
   c = 0.5*c +0.5*c_new;
   k_prime = 0.5*k_prime +0.5*k_prime_new;
end
end

