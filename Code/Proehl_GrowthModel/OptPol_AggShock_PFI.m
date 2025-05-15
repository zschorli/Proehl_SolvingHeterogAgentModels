% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file implements the policy function iteration for the model with
% aggregate shocks.
%__________________________________________________________________________
function [c,k_prime,diffs] ...
         = OptPol_AggShock_PFI(c,k_prime,basic_grid,StaticParams,Poly,grid_N,distrGrid,idx,wealth)

MargUt = @(c) c.^(-StaticParams.gamma);
MargUt_inv = @(c) c.^(-1./StaticParams.gamma);

diffs = 1e8.*ones(2,1);
iter=0;
agK_pr = zeros([size(k_prime,1),size(k_prime,2),size(k_prime,3)]);
c_pr = zeros([size(k_prime,1),4,4,size(k_prime,3)]);
c_new = zeros(size(c));
	
while diffs(max(iter,1)) > StaticParams.criter_k
   parfor m=1:size(k_prime,1)%
	   k_prime_large = interp1(StaticParams.kGrid_pol,squeeze(k_prime(m,:,:))',StaticParams.kGrid)';   
	   [agK_new,weights_new] = calcAggregates(k_prime_large,squeeze(basic_grid(m,:,:)),StaticParams,Poly);
	   agK_pr(m,:,:) = repmat(agK_new([1 1 2 2]),[1,length(StaticParams.kGrid_pol)]);
	   c_prime = @(z_idx,z,k,w)...
		   interp1(StaticParams.kGrid_pol,...
           getCApprox(grid_N,distrGrid,squeeze(c(:,2*z(1)+z(2)+1,:)),...
           idx,squeeze(w((z_idx>2)+1,z(1)+1,:))'),...
           max(StaticParams.k_min,k),'linear','extrap');
	   for z_idx=1:4
		   c_pr(m,z_idx,:,:) = [c_prime(z_idx,[0,0],squeeze(k_prime(m,z_idx,:))',weights_new);...
                                c_prime(z_idx,[0,1],squeeze(k_prime(m,z_idx,:))',weights_new);...
                                c_prime(z_idx,[1,0],squeeze(k_prime(m,z_idx,:))',weights_new);...
                                c_prime(z_idx,[1,1],squeeze(k_prime(m,z_idx,:))',weights_new)] ;
	   end
   end			
   rate_pr = @(z_idx) StaticParams.alpha.*...
				repmat([(1-StaticParams.delta_a)*ones(2,1);(1+StaticParams.delta_a)*ones(2,1)]',[size(k_prime,1),1,length(StaticParams.kGrid_pol)])...
			   .*(repmat(agK_pr(:,z_idx,:),[1,4,1])...
			   ./(StaticParams.l_bar.*repmat([StaticParams.er(1)*ones(2,1);StaticParams.er(2)*ones(2,1)]',[size(k_prime,1),1,length(StaticParams.kGrid_pol)])))...
			   .^(StaticParams.alpha-1);      
   for i=1:4
       c_new(:,i,:) = min(wealth(:,i,:),MargUt_inv(StaticParams.beta...
				    .*sum(repmat(StaticParams.P(i,:),[size(k_prime,1),1,length(StaticParams.kGrid_pol)])...
                    .*(1-StaticParams.delta+rate_pr(i))...
				    .*MargUt(squeeze(c_pr(:,i,:,:))),2)));
   end
   k_prime_new = max(StaticParams.k_min,wealth-c_new);
   iter = iter+1;
   diffs(iter) = max([max(max(squeeze(abs(k_prime_new(1,:,:)-k_prime(1,:,:))),[],2)),...
			          max(max(squeeze(abs(k_prime_new(2,:,:)-k_prime(2,:,:))),[],2))]);
   c = 0.5*c +0.5*c_new;
   k_prime = 0.5*k_prime +0.5*k_prime_new;
end
end

