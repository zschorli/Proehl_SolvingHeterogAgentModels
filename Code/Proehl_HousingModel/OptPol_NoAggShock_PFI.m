% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This file implements the policy function iteration for the model without
% aggregate shocks.
%__________________________________________________________________________
function [c,h,a_prime,diffs] ...
         = OptPol_NoAggShock_PFI(bool_gb,c,h,a_prime,StaticParams)

diffs = 1e8.*ones(2,1);
iter=0;
update = 0.1;
aux = zeros(size(c));
options = optimoptions('fsolve','Display','off');

innerAggr = @(c,h) c.^((StaticParams.eps-1)/StaticParams.eps)...
               +StaticParams.eta^(1/StaticParams.eps).*...
               h.^((StaticParams.eps-1)/StaticParams.eps);
Aggr = @(c,h) (innerAggr(c,h)).^(StaticParams.eps/(StaticParams.eps-1));
Aggr_derC = @(c,h) (innerAggr(c,h)).^(StaticParams.eps/(StaticParams.eps-1)-1)...
               .*c.^(-1/StaticParams.eps);
MargUt = @(g) g.^(-StaticParams.gamma);
MargUt_inv = @(g) g.^(-1./StaticParams.gamma);

if bool_gb
   P = StaticParams.P_g;
else
   P = StaticParams.P_b;
end   
wealth = StaticParams.wealth(bool_gb);
% OptC_mult results from inserting unconstrained optimal h into utility function
if 1/StaticParams.eps==StaticParams.gamma
    OptC_mult = ones(size(StaticParams.uq_large)); 
else
    OptC_mult = (1+StaticParams.eta.*StaticParams.uq_large.^(1-StaticParams.eps))...
                .^((1-StaticParams.eps*StaticParams.gamma)./(StaticParams.eps-1)); 
end
a_large = repmat(repmat(StaticParams.aGrid_pol,[2,1]),[1,1,StaticParams.Np,StaticParams.NA]);
normalized_a = a_large./StaticParams.q_large;
F_h = @(c) min((StaticParams.eta.*StaticParams.uq_large.^(-StaticParams.eps).*c)...
          ,(1/(1-StaticParams.theta).*normalized_a));
F_h_idx = @(c,idx) min((StaticParams.eta.*StaticParams.uq_large(idx).^(-StaticParams.eps).*c)...
          ,(1/(1-StaticParams.theta).*normalized_a(idx)));
c_max = c;
parfor i=1:numel(c)
    c_max(i) = fsolve(@(c) wealth(i)-StaticParams.uq_large(i).*F_h_idx(c,i)-c...
                     ,wealth(i),options);
end
	
while diffs(max(iter,1)) > StaticParams.criter_a
   MargUt_pr = @(z_ag) StaticParams.beta.*StaticParams.R_large(z_ag+1,:,:,:)...
				.*sum(repmat(P(z_ag+1,:)',[1,length(StaticParams.aGrid_pol),StaticParams.Np,StaticParams.NA])...
                .*MargUt(Aggr(max(eps,fct_a(StaticParams.aGrid_pol,c,a_prime([z_ag+1 z_ag+1],:,:,:))),...
                              max(eps,fct_a(StaticParams.aGrid_pol,h,a_prime([z_ag+1 z_ag+1],:,:,:)))))...
                .*Aggr_derC(max(eps,fct_a(StaticParams.aGrid_pol,c,a_prime([z_ag+1 z_ag+1],:,:,:))),...
                            max(eps,fct_a(StaticParams.aGrid_pol,h,a_prime([z_ag+1 z_ag+1],:,:,:)))),1);
   aux(1,:,:,:) = MargUt_inv(MargUt_pr(0)./OptC_mult(1,:,:,:));
   aux(2,:,:,:) = MargUt_inv(MargUt_pr(1)./OptC_mult(2,:,:,:));
   % c_new = aux may not be monotone as the constraint introduces 
   % oscillations. The following line corrects for this. This is crucial 
   % for convergence
   aux2 = cummin(aux,2,"reverse");
   c_new = max(0,min(c_max,aux2));
    
   iter = iter+1;
   diffs(iter,1) = max(max(max(max(abs(c_new-c)))));
   diffs(iter,2) = max(max(max(max(abs(c_new-max(0,min(c_max,aux)))))));
   c = (1-update)*c + update*c_new;
   h = F_h(c);
   a_prime = max(0,wealth-c-StaticParams.uq_large.*h);
   a_prime(a_prime<10*eps) = 0; 
end
end