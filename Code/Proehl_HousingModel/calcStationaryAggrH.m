% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This file computes the stationary aggregate labor for given policy
% functions
%__________________________________________________________________________
function H = calcStationaryAggrH(a_prime_approx,h_approx,P,p,StaticParams)
[~,H] = calcStationaryDistr(a_prime_approx,h_approx,P,p,StaticParams);
end