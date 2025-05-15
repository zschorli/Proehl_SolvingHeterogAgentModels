% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file updates labor supply in the policy function iteration in line
% with the FOCs
%__________________________________________________________________________
function [f,J] = solveLaborEqu(l,aux_unconstr,StaticParams)
   % log(aux_constr) = -gamma*log(1-l)+(rho*(1-(1-gamma)*theta))*log(l)
   f = log(aux_unconstr)+StaticParams.gamma.*log(1-l)-(StaticParams.rho.*(1-(1-StaticParams.gamma).*StaticParams.theta)).*log(l);
   j_aux = -StaticParams.gamma./(1-l)-(StaticParams.rho.*(1-(1-StaticParams.gamma).*StaticParams.theta))./(l);
   J = sparse(1:numel(j_aux),1:numel(j_aux),j_aux);
end