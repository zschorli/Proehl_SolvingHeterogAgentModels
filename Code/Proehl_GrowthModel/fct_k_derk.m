% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function interpolates along the individual capital.
%__________________________________________________________________________

function der_k = fct_k_derk(kGrid_pol,k_prime,grid) 

der_k = zeros(size(k_prime));
h = 1e-6;
for i=1:size(k_prime,3)
for j=1:size(k_prime,1)
	der_k(j,:,i) = (interp1(kGrid_pol,squeeze(k_prime(j,:,i)),squeeze(grid(j,:,i))+h,'linear','extrap')...
                   -interp1(kGrid_pol,squeeze(k_prime(j,:,i)),squeeze(grid(j,:,i))-h,'linear','extrap'))./(2*h);
end
end