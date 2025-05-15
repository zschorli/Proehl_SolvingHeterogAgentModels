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
function [c,l,k_prime,diffs] ...
         = OptPol_NoAggShock_PFI(bool_gb,c,l,k_prime,StaticParams)
 
diffs = 1e8.*ones(1,2);
iter=0;
update = 0.5;
options = optimoptions(@fsolve,'Display','off','SpecifyObjectiveGradient',true,...
          'MaxIterations',1000,'FunctionTolerance',StaticParams.criter_k,'StepTolerance',1e-10);

c_new = zeros(size(c));
l_new = zeros(size(l));
if bool_gb
   P = StaticParams.P_g;
else
   P = StaticParams.P_b;
end   
MargUt = @(u) u.^(-StaticParams.gamma);
rate = repmat(StaticParams.rate(bool_gb,StaticParams.KGrid,StaticParams.LGrid),[1,numel(StaticParams.kGrid_pol),1,1,numel(StaticParams.LtaxGrid)]);
% u = c^theta*(1-l)^(1-theta) = (1-l)*l^(-rho*theta)*Q^theta if the FOC for l holds
% Note that Q = l^(rho)*c/(1-l) is independent of l and c
Q_unconstr = repmat(StaticParams.Q(bool_gb,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid),[1,numel(StaticParams.kGrid_pol),1,1,1]);


while diffs(max(iter,1)) > StaticParams.criter_k
   l_new(1,:,:,:,:) = StaticParams.lu;
   c_pr = max(eps,fct_k(StaticParams.kGrid_pol,c,k_prime([1 1],:,:,:,:)));
   l_pr = max(0,min(1-eps,fct_k(StaticParams.kGrid_pol,l,k_prime([1 1],:,:,:,:))));
   l_pr(2,:,:,:,:) = max(StaticParams.l_min,l_pr(2,:,:,:,:));
   l_pr(1,:,:,:,:) = StaticParams.lu;
   RHS_FOCc = StaticParams.beta...
			  .*sum(repmat(P(1,:)',[1,length(StaticParams.kGrid_pol),length(StaticParams.KGrid),length(StaticParams.LGrid),length(StaticParams.LtaxGrid)])...
                    .*(1-StaticParams.delta+rate(1,:,:,:,:))...
			        .*MargUt(c_pr.^StaticParams.theta.*(1-l_pr).^(1-StaticParams.theta))...
                    .*c_pr.^(StaticParams.theta-1).*(1-l_pr).^(1-StaticParams.theta),1);
   c_new(1,:,:,:,:) = (RHS_FOCc.*(1-l_new(1,:,:,:,:)).^(-(1-StaticParams.theta).*(1-StaticParams.gamma)))...
                      .^(1/((1-StaticParams.gamma).*StaticParams.theta-1));
   
   c_pr = max(eps,fct_k(StaticParams.kGrid_pol,c,k_prime([2 2],:,:,:,:)));
   l_pr = min(1-eps,fct_k(StaticParams.kGrid_pol,l,k_prime([2 2],:,:,:,:)));
   l_pr(2,:,:,:,:) = max(StaticParams.l_min,l_pr(2,:,:,:,:));
   l_pr(1,:,:,:,:) = StaticParams.lu;
   RHS_FOCc = StaticParams.beta...
		      .*sum(repmat(P(2,:)',[1,length(StaticParams.kGrid_pol),length(StaticParams.KGrid),length(StaticParams.LGrid),length(StaticParams.LtaxGrid)])...
		            .*(1-StaticParams.delta+rate(2,:,:,:,:))...
                    .*MargUt(c_pr.^StaticParams.theta.*(1-l_pr).^(1-StaticParams.theta))...
                    .*c_pr.^(StaticParams.theta-1).*(1-l_pr).^(1-StaticParams.theta),1);
   aux_unconstr = RHS_FOCc.*Q_unconstr(2,:,:,:,:).^(1-(1-StaticParams.gamma)*StaticParams.theta);

   % u = (1-l)*l^(-rho*theta)*Q^theta such that
   % aux_unconstr = MargUt(1-l)*l^(rho*(1-(1-gamma)*theta))
   % fsolve is faster solving this equation but suffers from instabilities
   % in the constrained region when taxes are progressive, I use fzero instead
   a_match = reshape(aux_unconstr,[],1);
   aux_l = zeros(numel(aux_unconstr),1);
   a_start = reshape(squeeze(l_new(2,:,:,:,:)),[],1);
   if StaticParams.rho==0
       [aux_l,~,flag] = fsolve(@(l) solveLaborEqu(l,a_match,StaticParams),max(eps,a_start),options);
   elseif iter>=5
       % identify unconstrained region
       a_bool = logical(1-(a_start==StaticParams.l_min));
       a_start = max(eps,a_start-0.1);
       % a_bool(:) = true;
       [aux_l(a_bool),~,flag] = fsolve(@(l) solveLaborEqu(l,a_match(a_bool),StaticParams),a_start(a_bool),options);
       if flag<0
          parfor i=1:numel(aux_l)
              aux_l(i) = fzero(@(l) solveLaborEqu(l,a_match(i),StaticParams),[eps,1-eps]);
          end
       elseif sum(1-a_bool)>0
          idx = 1:numel(aux_unconstr);
          idx = idx(logical(1-a_bool));
          aux_sol = zeros(size(idx));
          aux_match = a_match(idx);
          parfor i=1:sum(1-a_bool)
            aux_sol(i) = fzero(@(l) solveLaborEqu(l,aux_match(i),StaticParams),[eps,1-eps]);
          end 
          aux_l(idx) = aux_sol;
       end
   else
       parfor i=1:numel(aux_l)
           aux_l(i) = fzero(@(l) solveLaborEqu(l,a_match(i),StaticParams),[eps,1-eps]);
       end
   end

   aux_l = max(StaticParams.l_min,min(1,reshape(real(aux_l),size(aux_unconstr))));
   l_new(2,:,:,:,:) = aux_l;
   % c = (1-l)*l^(-rho)*Q when l>l_min, 
   % otherwise c=(RHS_FOCc*(1-l_min)^(-(1-theta)*(1-gamma)))^(1/((1-gamma)*theta-1)) when l=l_min
   aux_c = (RHS_FOCc.*(1-StaticParams.l_min).^(-(1-StaticParams.theta).*(1-StaticParams.gamma)))...
           .^(1/((1-StaticParams.gamma).*StaticParams.theta-1));
   aux_Q = Q_unconstr(2,:,:,:,:);
   aux_c(aux_l>StaticParams.l_min) = (1-aux_l(aux_l>StaticParams.l_min)).*aux_l(aux_l>StaticParams.l_min).^(-StaticParams.rho).*aux_Q(aux_l>StaticParams.l_min);
   c_new(2,:,:,:,:) = aux_c;

   wealth = StaticParams.wealth(bool_gb,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,StaticParams.kGrid_pol,l_new);
   k_prime_new = max(StaticParams.k_min,wealth-c_new);
   c_new = min(wealth-StaticParams.k_min,c_new);

   iter = iter+1;
   diffs(iter,:) = [mean(mean(mean(mean(mean(abs(c_new-c)))))),mean(mean(mean(mean(mean(abs(l_new-l))))))];
   
   c = (1-update)*c +update*c_new;
   l = (1-update)*l +update*l_new;
   k_prime = (1-update)*k_prime +update*k_prime_new;
end
end