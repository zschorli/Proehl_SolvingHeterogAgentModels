% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This file implements the aggregate asset demand/supply to find the 
% equilibrium prices.
%__________________________________________________________________________
function [A,pdf_cond] = AggrAssetAlloc(p,k_prime,pdf_cond,StaticParams,KGrid,bool_gb,bool_stationary,Poly)
[KGrid2,idx] = unique(KGrid);
if bool_stationary
    k_prime_approx = @(z_idx,k) interpn(StaticParams.kGrid_pol,KGrid2,squeeze(k_prime(z_idx,:,idx)),k,p);
    iter = 1;
    diffs = 1e8*ones(2,1);
    while (diffs(end)) > 5e-7 && iter<5000
        cdf_cond_new = cummax(min(1,cumsum(pdf_cond,2)),2);
        cdf_cond_new(1,cdf_cond_new(1,:)==max(cdf_cond_new(1,:))) = 1;
        cdf_cond_new(2,cdf_cond_new(2,:)==max(cdf_cond_new(2,:))) = 1;
        k_prime = sort([k_prime_approx(1,StaticParams.kGrid_pol)';...
                        k_prime_approx(2,StaticParams.kGrid_pol)'],2);
        pdf_cond_new = calcEoPDistr(max(StaticParams.k_min,interp1(StaticParams.kGrid_pol,k_prime',StaticParams.kGrid)'),...
                       cdf_cond_new,StaticParams.kGrid);
        pdf_cond_new = (repmat(1./(StaticParams.p(2*bool_gb+(1:2))'*StaticParams.P(2*bool_gb+(1:2),2*bool_gb+(1:2)))',[1,2])...
                       .*StaticParams.P(2*bool_gb+(1:2),2*bool_gb+(1:2))')...
                       *(repmat(StaticParams.p(2*bool_gb+(1:2)),[1,size(pdf_cond_new,2)])...
                       .*pdf_cond_new);
        if bool_stationary
            diffs(iter) = norm(pdf_cond_new-pdf_cond,2);
        else
            diffs(iter) = 0;
        end
        pdf_cond_new = pdf_cond_new.*(pdf_cond_new>eps);
        pdf_cond = pdf_cond_new./repmat(sum(pdf_cond_new,2),[1,size(pdf_cond_new,2)]);
        iter = iter+1;
    end
    A = 2*StaticParams.p(2*bool_gb+(1:2))'*pdf_cond*StaticParams.kGrid';
else
    if size(k_prime,1)<4
        k_prime = repmat(k_prime,[2,1,1]);
    end
    if length(KGrid2)>1
        k_prime = interpn(1:2,StaticParams.kGrid_pol,KGrid2,squeeze(k_prime(2*bool_gb+(1:2),:,idx)),1:2,StaticParams.kGrid,p);
    else
        k_prime = interp1(StaticParams.kGrid_pol,squeeze(k_prime(2*bool_gb+(1:2),:))',StaticParams.kGrid)';
    end
    A = calcAggregates(k_prime,pdf_cond,StaticParams,Poly);
    A = A(bool_gb+1);
end
end