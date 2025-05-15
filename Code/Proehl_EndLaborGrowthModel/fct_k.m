% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function interpolates along the individual capital.
%__________________________________________________________________________

function k = fct_k(kGrid_pol,k_prime,grid) 

k = zeros(size(k_prime));

for m=1:size(k_prime,5)
for l=1:size(k_prime,4)
for i=1:size(k_prime,3)
for j=1:size(k_prime,1)
	k(j,:,i,l,m) = interp1(kGrid_pol,squeeze(k_prime(j,:,i,l,m)),squeeze(grid(j,:,i,l,m)),'linear','extrap');
end
end
end
end