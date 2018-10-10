% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This function implements the inverse of the c.d.f. needed to compute the
% weights of the polynomial chaos expansion.
%__________________________________________________________________________

function kGrid_new = projectDistr(kGrid,cdf,cdf_new)
if isempty(find(cdf==0,1,'last'))==0
   count = find(cdf==0,1,'last');
   kGrid = kGrid(count:end);
   cdf = cdf(count:end);
end
[cdf,idx] = unique(cdf,'first','legacy');
kGrid = kGrid(idx);
if length(cdf)==1
    cdf=[0,cdf];
    kGrid=[kGrid-1e-5,kGrid];
end
kGrid_new = interp1(cdf,kGrid,max(cdf(1),cdf_new));
end