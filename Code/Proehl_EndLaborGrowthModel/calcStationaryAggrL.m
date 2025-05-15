% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file computes the stationary aggregate labor for given policy
% functions
%__________________________________________________________________________
function L = calcStationaryAggrL(k_prime_approx,l_approx,P,p,StaticParams)
[~,L] = calcStationaryDistr(k_prime_approx,l_approx,P,p,StaticParams);
end