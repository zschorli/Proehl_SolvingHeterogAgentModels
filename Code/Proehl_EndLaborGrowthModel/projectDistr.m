% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pr�hl
%
% AUTHOR Elisabeth Pr�hl, University of Amsterdam
% DATE May 2025
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
kGrid_new = interp1(cdf,kGrid,max(cdf(1),cdf_new));
end