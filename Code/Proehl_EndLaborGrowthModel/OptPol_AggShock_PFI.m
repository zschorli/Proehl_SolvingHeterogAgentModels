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
function [c,l,k_prime,agL,agLtax,diffs] ...
         = OptPol_AggShock_PFI(c,l,k_prime,agL,agLtax,agK,basic_grid,StaticParams,Poly,grid_N,distrGrid,distrIdx)

MargUt = @(c) c.^(-StaticParams.gamma);
% Note that Q = l^(rho)*c/(1-l) is independent of l and c
Q_unconstr = zeros(size(k_prime));
for m=1:size(k_prime,1)
    Q_unconstr(m,1:2,:) = repmat(StaticParams.Q(0,agK(m),agL(m,1),agL(m,1)./agLtax(m,1)),[1,size(k_prime,3)]);
    Q_unconstr(m,3:4,:) = repmat(StaticParams.Q(1,agK(m),agL(m,2),agL(m,2)./agLtax(m,2)),[1,size(k_prime,3)]);
end
wealth = zeros(size(k_prime));

diffs = 1e8.*ones(1,3);
iter=0;
options = optimoptions(@fsolve,'Display','off','SpecifyObjectiveGradient',true,...
          'MaxIterations',1000,'FunctionTolerance',StaticParams.criter_k,'StepTolerance',1e-10);
update = 0.3;

agK_pr = zeros([size(k_prime,1),4,4,size(k_prime,3)]);
agL_pr = zeros([size(k_prime,1),4,4,size(k_prime,3)]);
c_pr = zeros([size(k_prime,1),4,4,size(k_prime,3)]);
l_pr = zeros([size(k_prime,1),4,4,size(k_prime,3)]);
c_new = zeros(size(c));
l_new = zeros(size(l));
agL_new = zeros(size(agL));
agLtax_new = zeros(size(agL));

while (max(diffs(max(iter,1),:)) > StaticParams.criter_k) && (iter<=2000)
   parfor m=1:size(k_prime,1)%
	   k_prime_large = max(StaticParams.k_min,interp1(StaticParams.kGrid_pol,squeeze(k_prime(m,:,:))',StaticParams.kGrid)');   
	   [agK_new,weights_new] = calcAggregates(k_prime_large,squeeze(basic_grid(m,:,:)),StaticParams,Poly);
	   agK_pr(m,:,:,:) = repmat(agK_new([1 1 2 2]),[1,4,length(StaticParams.kGrid_pol)]);
	   l_prime = @(z_idx,z,k,w)...
		   interp1(StaticParams.kGrid_pol,...
           getCApprox(grid_N,distrGrid,squeeze(l(:,2*z(1)+z(2)+1,:)),...
           distrIdx,squeeze(w((z_idx>2)+1,z(1)+1,:))'),...
           max(StaticParams.k_min,k),'linear','extrap');
	   c_prime = @(z_idx,z,k,w)...
		   interp1(StaticParams.kGrid_pol,...
           getCApprox(grid_N,distrGrid,squeeze(c(:,2*z(1)+z(2)+1,:)),...
           distrIdx,squeeze(w((z_idx>2)+1,z(1)+1,:))'),...
           max(StaticParams.k_min,k),'linear','extrap');
	   for z_idx=1:4
		   l_pr_large = max(StaticParams.l_min,min(1,[l_prime(z_idx,[0,0],squeeze(k_prime_large(z_idx,:)),weights_new);...
                               l_prime(z_idx,[0,1],squeeze(k_prime_large(z_idx,:)),weights_new);...
                               l_prime(z_idx,[1,0],squeeze(k_prime_large(z_idx,:)),weights_new);...
                               l_prime(z_idx,[1,1],squeeze(k_prime_large(z_idx,:)),weights_new)]));
           [~,~,~,agL_pr_new] = calcAggregates(l_pr_large,squeeze(basic_grid(m,:,:)),StaticParams,Poly);
           agL_pr(m,z_idx,:,:) = repmat(agL_pr_new([1 1 2 2]),[1,length(StaticParams.kGrid_pol)]);
		   c_pr(m,z_idx,:,:) = max(eps,[c_prime(z_idx,[0,0],squeeze(k_prime(m,z_idx,:))',weights_new);...
                                c_prime(z_idx,[0,1],squeeze(k_prime(m,z_idx,:))',weights_new);...
                                c_prime(z_idx,[1,0],squeeze(k_prime(m,z_idx,:))',weights_new);...
                                c_prime(z_idx,[1,1],squeeze(k_prime(m,z_idx,:))',weights_new)]);
           if (z_idx==2) || (z_idx==4)
		       l_pr(m,z_idx,:,:) = [l_prime(z_idx,[0,0],squeeze(k_prime(m,z_idx,:))',weights_new);...
                                    l_prime(z_idx,[0,1],squeeze(k_prime(m,z_idx,:))',weights_new);...
                                    l_prime(z_idx,[1,0],squeeze(k_prime(m,z_idx,:))',weights_new);...
                                    l_prime(z_idx,[1,1],squeeze(k_prime(m,z_idx,:))',weights_new)];
               l_pr(m,z_idx,:,:) = min(1,max(StaticParams.l_min,l_pr(m,z_idx,:,:)));
           else
               l_pr(m,z_idx,:,:) = StaticParams.lu;
           end
	   end
   end
   rate_pr = @(z_idx) StaticParams.alpha.*...
				repmat([(1-StaticParams.delta_a)*ones(2,1);(1+StaticParams.delta_a)*ones(2,1)]',[size(k_prime,1),1,length(StaticParams.kGrid_pol)])...
			   .*(squeeze(agK_pr(:,z_idx,:,:))...
			   ./(StaticParams.l_bar.*repmat([StaticParams.er(1)*ones(2,1);StaticParams.er(2)*ones(2,1)]',[size(k_prime,1),1,length(StaticParams.kGrid_pol)])...
                 .*squeeze(agL_pr(:,z_idx,:,:))))...
			   .^(StaticParams.alpha-1); 
   for i=[2,4]
       RHS_FOCc = StaticParams.beta...
  			      .*sum(repmat(StaticParams.P(i,:),[size(k_prime,1),1,length(StaticParams.kGrid_pol)])...
                        .*(1-StaticParams.delta+rate_pr(i))...
				        .*MargUt(squeeze(c_pr(:,i,:,:)).^StaticParams.theta.*squeeze(1-l_pr(:,i,:,:)).^(1-StaticParams.theta))...
                        .*squeeze(c_pr(:,i,:,:)).^(StaticParams.theta-1).*squeeze(1-l_pr(:,i,:,:)).^(1-StaticParams.theta),2);
       aux_unconstr = RHS_FOCc.*Q_unconstr(:,i,:).^(1-(1-StaticParams.gamma)*StaticParams.theta);
       % u = (1-l)*l^(-rho*theta)*Q^theta such that
       % aux_unconstr = MargUt(1-l)*l^(rho*(1-(1-gamma)*theta))
       a_match = reshape(aux_unconstr,[],1);
       a_start = reshape(l(:,i,:),[],1);
       aux_l = zeros(numel(aux_unconstr),1);
       if StaticParams.rho==0
           idx = 1:length(a_start);
           split=floor(length(idx)/6)+1;
           for j=1:6
                idx_select = idx((j-1)*split+1:min(length(idx),j*split));
                [aux_l(idx_select),~,flag] = fsolve(@(l) solveLaborEqu(l,a_match(idx_select),StaticParams),a_start(idx_select),options);
           end
       else 
           % identify unconstrained region
           a_bool = logical(1-(a_start==StaticParams.l_min));
           a_start = max(eps,a_start-0.1);
           idx = 1:length(a_bool);
           idx = idx(a_bool);
           split=floor(length(idx)/6)+1;
           flag = 1;
           for j=1:6
                idx_select = idx((j-1)*split+1:min(length(idx),j*split));
                [aux_l(idx_select),~,aux_flag] = fsolve(@(l) solveLaborEqu(l,a_match(idx_select),StaticParams),a_start(idx_select),options);
                flag = min(flag,aux_flag);
           end
           if flag<0
              parfor j=1:numel(aux_l)
                  aux_l(j) = fzero(@(l) solveLaborEqu(l,a_match(j),StaticParams),[eps,1-eps]);
              end
           elseif sum(1-a_bool)>0
              idx = 1:numel(aux_unconstr);
              idx = idx(logical(1-a_bool));
              aux_sol = zeros(size(idx));
              aux_match = a_match(idx);
              parfor j=1:sum(1-a_bool)
                aux_sol(j) = fzero(@(l) solveLaborEqu(l,aux_match(j),StaticParams),[eps,1-eps]);
              end 
              aux_l(idx) = aux_sol;
           end
       end
       aux_l = max(StaticParams.l_min,min(1,reshape(real(aux_l),size(aux_unconstr))));
       l_new(:,i,:) = aux_l;
       % c = (1-l)*l^(-rho)*Q when l>l_min, otherwise
       % c=(RHS_FOCc*(1-l_min)^(-(1-theta)*(1-gamma)))^(1/((1-gamma)*theta-1)) when l=l_min
       aux_c = (RHS_FOCc.*(1-StaticParams.l_min).^(-(1-StaticParams.theta).*(1-StaticParams.gamma)))...
               .^(1/((1-StaticParams.gamma).*StaticParams.theta-1));
       aux_Q = Q_unconstr(:,i,:);
       aux_c(aux_l>StaticParams.l_min) = (1-aux_l(aux_l>StaticParams.l_min)).*aux_l(aux_l>StaticParams.l_min).^(-StaticParams.rho).*aux_Q(aux_l>StaticParams.l_min);
       c_new(:,i,:) = aux_c;
   end
   parfor m=1:size(k_prime,1)
        l_large = interp1(StaticParams.kGrid_pol,squeeze(l_new(m,:,:))',StaticParams.kGrid)';
        [~,~,~,agL_new(m,:)] = calcAggregates(l_large,squeeze(basic_grid(m,:,:)),StaticParams,Poly);
        [~,~,~,agLtax_new(m,:)] = calcAggregates(l_large.^(1-StaticParams.rho),squeeze(basic_grid(m,:,:)),StaticParams,Poly);
   end
   l_new(:,[1,3],:) = StaticParams.lu;
   for i=[1,3]
       RHS_FOCc = StaticParams.beta...
  			      .*sum(repmat(StaticParams.P(i,:),[size(k_prime,1),1,length(StaticParams.kGrid_pol)])...
                        .*(1-StaticParams.delta+rate_pr(i))...
				        .*MargUt(squeeze(c_pr(:,i,:,:)).^StaticParams.theta.*squeeze(1-l_pr(:,i,:,:)).^(1-StaticParams.theta))...
                        .*squeeze(c_pr(:,i,:,:)).^(StaticParams.theta-1).*squeeze(1-l_pr(:,i,:,:)).^(1-StaticParams.theta),2);
       c_new(:,i,:) = (RHS_FOCc...
                    .*(1-l_new(:,i,:)).^(-(1-StaticParams.theta).*(1-StaticParams.gamma)))...
                    .^(1/((1-StaticParams.gamma).*StaticParams.theta-1));
   end
   parfor m=1:size(k_prime,1)%
        wealth(m,:,:) = [StaticParams.wealth(0,agK(m),agL_new(m,1),agL_new(m,1)./agLtax_new(m,1),StaticParams.kGrid_pol,squeeze(l_new(m,1:2,:)));...
                         StaticParams.wealth(1,agK(m),agL_new(m,2),agL_new(m,2)./agLtax_new(m,2),StaticParams.kGrid_pol,squeeze(l_new(m,3:4,:)))];
   end
   k_prime_new = max(StaticParams.k_min,wealth-c_new);
   c_new = min(wealth-StaticParams.k_min,c_new);

   iter = iter+1;
   diffs(iter,:) = [mean(mean(mean(abs(c_new-c)))),mean(mean(mean(abs(l_new-l)))),mean(mean(abs(agL_new-agL)))];

   c = (1-update)*c +update*c_new;
   l = (1-update)*l +update*l_new;
   agL = (1-update)*agL +update*agL_new;
   agLtax = (1-update)*agLtax +update*agLtax_new;
   k_prime = (1-update)*k_prime +update*k_prime_new;
end
end

