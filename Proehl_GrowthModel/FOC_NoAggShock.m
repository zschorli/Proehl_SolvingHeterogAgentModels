% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This function computes the first-order conditions to solve for in the
% proximal point algorithm without aggregate shocks.
%__________________________________________________________________________
function [f,J] = FOC_NoAggShock(sol,c,k_prime,y,i,bool_gb,m,StaticParams,preSol)

h = 1e-6;
sol_k = reshape(sol,size(k_prime'))';


c_prime = @(z,k) interp1(StaticParams.kGrid_pol,c(z+1,:),max(StaticParams.k_min,k),'linear','extrap');  
c_pr = zeros([2,size(c)]);
c_pr_der = zeros([2,size(c)]);
for z_idx=1:2
    c_pr(z_idx,:,:)=[c_prime(0,sol_k(z_idx,:));...
                     c_prime(1,sol_k(z_idx,:))] ;
    c_pr_der(z_idx,:,:)=[(c_prime(0,sol_k(z_idx,:)+h)-c_prime(0,sol_k(z_idx,:)-h))./(2*h);...
                         (c_prime(1,sol_k(z_idx,:)+h)-c_prime(1,sol_k(z_idx,:)-h))./(2*h)] ;
end        

f_foc_k = @(z_idx)...
       (squeeze(preSol.wealth(bool_gb*2+z_idx,:,i))-sol_k(z_idx,:)).^(-StaticParams.gamma)...
       -StaticParams.beta*(...
        (StaticParams.P(z_idx+2*bool_gb,2*bool_gb+(1:2))./sum(StaticParams.P(z_idx+2*bool_gb,2*bool_gb+(1:2))))*...
        ((1-StaticParams.delta+repmat(squeeze(preSol.rate(2*bool_gb+z_idx,:,i)),[2,1]))...
        .*squeeze(c_pr(z_idx,:,:)).^(-StaticParams.gamma)))...
       +1/m*(sol_k(z_idx,:)-k_prime(z_idx,:))...
       +(sol_k(z_idx,:)<=y(z_idx,:)/m).*(-y(z_idx,:)+m*(sol_k(z_idx,:)));
f_foc_k_dk = @(z_idx)...
       (StaticParams.gamma).*(squeeze(preSol.wealth(bool_gb*2+z_idx,:,i))-sol_k(z_idx,:)).^(-StaticParams.gamma-1)...
       -StaticParams.beta*(...
        (StaticParams.P(z_idx+2*bool_gb,2*bool_gb+(1:2))/sum(StaticParams.P(z_idx+2*bool_gb,2*bool_gb+(1:2))))...
        *((1-StaticParams.delta+repmat(squeeze(preSol.rate(2*bool_gb+z_idx,:,i)),[2,1]))...
        .*(-StaticParams.gamma).*squeeze(c_pr(z_idx,:,:)).^(-StaticParams.gamma-1)...
        .*squeeze(c_pr_der(z_idx,:,:))))...
       +1/m+(sol_k(z_idx,:)<=y(z_idx,:)/m).*m;
   
foc_k = [f_foc_k(1);f_foc_k(2)].*squeeze(preSol.wealth((1:2)+2*bool_gb,:,i));
foc_k_dk = [f_foc_k_dk(1);f_foc_k_dk(2)].*squeeze(preSol.wealth((1:2)+2*bool_gb,:,i));

% function value of the Euler equation
f = reshape(foc_k',1,numel(foc_k))';

% Jacobian of the Euler equation
if nargout>1    
    J = sparse(1:length(sol),1:length(sol),reshape(foc_k_dk',numel(foc_k_dk),1));
end
end

