% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function computes the first-order conditions to solve for in the
% proximal point algorithm without aggregate shocks.
%__________________________________________________________________________
function [f,J] = FOC_NoAggShock(sol,k_prime,y,i,bool_gb,m,StaticParams)

sol_k = reshape(sol,size(k_prime));
if bool_gb
   P = StaticParams.P_g;
else
   P = StaticParams.P_b;
end 
ret = @(k) 1-StaticParams.delta+repmat(StaticParams.rate(bool_gb,StaticParams.KGrid(i)),[1,numel(k),1]);
c_pr = @(k) StaticParams.wealth(bool_gb,StaticParams.KGrid(i),k)...
            -fct_k(StaticParams.kGrid_pol,sol_k,repmat(k,[2,1]));
c_pr_der = @(k) ret(k)-fct_k_derk(StaticParams.kGrid_pol,sol_k,repmat(k,[2,1]));
MargUt = @(c) c.^(-StaticParams.gamma);
MargUt_inv = @(c) c.^(-1./StaticParams.gamma);
MargUt_der = @(c) (-StaticParams.gamma).*c.^(-StaticParams.gamma-1);
MargUt_inv_der = @(c) (-1./StaticParams.gamma).*c.^(-1./StaticParams.gamma-1);

foc_k = -(c_pr(StaticParams.kGrid_pol))...
       +MargUt_inv(StaticParams.beta.*ret(StaticParams.kGrid_pol)...
        .*[P(1,:)*MargUt(c_pr(sol_k(1,:)));P(2,:)*MargUt(c_pr(sol_k(2,:)))])...
       +1./m.*(sol_k-k_prime)...
       +(sol_k<=y./m).*(-y-m.*sol_k)+(sol_k>y./m).*(-1./m.*y);
foc_k_dk = +1 ...
       +MargUt_inv_der(StaticParams.beta.*ret(StaticParams.kGrid_pol)...
        .*[P(1,:)*MargUt(c_pr(sol_k(1,:)));P(2,:)*MargUt(c_pr(sol_k(2,:)))])...
       .*(StaticParams.beta.*ret(StaticParams.kGrid_pol)...
        .*[P(1,:)*(MargUt_der(c_pr(sol_k(1,:))).*c_pr_der(sol_k(1,:)));...
           P(2,:)*(MargUt_der(c_pr(sol_k(2,:))).*c_pr_der(sol_k(2,:)))])...
       +1./m+(sol_k<=y./m).*(-m);

% function value of the Euler equation
f = reshape(foc_k,[],1);

% Jacobian of the Euler equation
if nargout>1    
    J = sparse(1:length(sol),1:length(sol),reshape(foc_k_dk,[],1));
end
end

