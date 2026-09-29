% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This function implements the inverse of the c.d.f. needed to compute the
% weights of the polynomial chaos expansion.
%__________________________________________________________________________

function aGrid_new = projectDistr(aGrid,cdf,cdf_new)
if isempty(find(cdf==0,1,'last'))==0
   count = find(cdf==0,1,'last');
   aGrid = aGrid(count:end);
   cdf = cdf(count:end);
end
[cdf,idx] = unique(cdf,'first','legacy');
aGrid = aGrid(idx);
if length(cdf)==1
    cdf=[0,cdf];
    aGrid=[aGrid-1e-5,aGrid];
end
aGrid_new = interp1(cdf,aGrid,max(cdf(1),cdf_new));
end