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

for i=1:size(k_prime,3)
for j=1:size(k_prime,1)
	k(j,:,i) = interp1(kGrid_pol,squeeze(k_prime(j,:,i)),squeeze(grid(j,:,i)),'linear','extrap');
end
end