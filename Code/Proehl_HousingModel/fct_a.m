% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This function interpolates along the individual net worth.
%__________________________________________________________________________

function a = fct_a(aGrid_pol,a_prime,grid) 

a = zeros(size(a_prime));

for m=1:size(a_prime,5)
for l=1:size(a_prime,4)
for i=1:size(a_prime,3)
for j=1:size(a_prime,1)
	a(j,:,i,l,m) = interp1(aGrid_pol,squeeze(a_prime(j,:,i,l,m)),squeeze(grid(j,:,i,l,m)),'linear','extrap');
end
end
end
end