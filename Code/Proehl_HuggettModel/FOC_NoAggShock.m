% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This function computes the first-order conditions to solve for in the
% proximal point algorithm without aggregate shocks.
%__________________________________________________________________________
function [f,J] = FOC_NoAggShock(sol,c,k_prime,y,y2,i,bool_gb,m,StaticParams,preSol)

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
if bool_gb
   P = StaticParams.P_g;
else
   P = StaticParams.P_b;
end   

f_foc_k = @(z_idx)...
       preSol.KGrid(i).*max(eps,squeeze(preSol.wealth(bool_gb*2+z_idx,:,i))-sol_k(z_idx,:).*preSol.KGrid(i)).^(-StaticParams.gamma)...
       -StaticParams.beta*(...
        P(z_idx,:)*...
        (squeeze(c_pr(z_idx,:,:)).^(-StaticParams.gamma)))...
       +1/m*(sol_k(z_idx,:)-k_prime(z_idx,:))...
       +((sol_k(z_idx,:)-StaticParams.k_min)<=y(z_idx,:)/m).*(-y(z_idx,:)+m*(sol_k(z_idx,:)-StaticParams.k_min))...
       +((squeeze(preSol.wealth(bool_gb*2+z_idx,:,i))-sol_k(z_idx,:).*preSol.KGrid(i))<=y2(z_idx,:)/m).*(-y2(z_idx,:)+m*(squeeze(preSol.wealth(bool_gb*2+z_idx,:,i))-sol_k(z_idx,:).*preSol.KGrid(i)));%...
   
foc_k = [f_foc_k(1);f_foc_k(2)].*squeeze(preSol.wealth((1:2)+2*bool_gb,:,i));

% function value of the Euler equation
f = reshape(foc_k',1,numel(foc_k))';

% Jacobian of the Euler equation
if nargout>1    
    f_foc_k_dk = @(z_idx)...
           (StaticParams.gamma).*preSol.KGrid(i)^2 ...
           .*max(eps,squeeze(preSol.wealth(bool_gb*2+z_idx,:,i))-sol_k(z_idx,:).*preSol.KGrid(i)).^(-StaticParams.gamma-1)...
           -StaticParams.beta*(...
            P(z_idx,:)...
            *((-StaticParams.gamma).*squeeze(c_pr(z_idx,:,:)).^(-StaticParams.gamma-1)...
            .*squeeze(c_pr_der(z_idx,:,:))))...
           +1/m+((sol_k(z_idx,:)-StaticParams.k_min)<=y(z_idx,:)/m).*m...
           -((squeeze(preSol.wealth(bool_gb*2+z_idx,:,i))-sol_k(z_idx,:).*preSol.KGrid(i))<=y2(z_idx,:)/m).*m*preSol.KGrid(i);
    foc_k_dk = [f_foc_k_dk(1);f_foc_k_dk(2)].*squeeze(preSol.wealth((1:2)+2*bool_gb,:,i));
    J = sparse(1:length(sol),1:length(sol),reshape(foc_k_dk',numel(foc_k_dk),1));
end
end

