% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This function calculates the law of motion from the beginning-of-period
% distribution to the end-of-period distribution after the agents apply 
% their capital policy.
%__________________________________________________________________________

function pdf_cond_new = calcEoPDistr(k_prime,cdf_cond_BoP,kGrid_pdf)

cdf_cond_new = zeros(size(kGrid_pdf));
pdf_cond_new = zeros(size(kGrid_pdf));

for j=1:size(k_prime,1)
    [kGrid_pr,idx] = unique(k_prime(j,:),'last','legacy');
    cdf = cdf_cond_BoP(j,idx);
    if kGrid_pr(1)>min(kGrid_pdf) && kGrid_pr(1)-1e6*eps>min(kGrid_pdf) 
        kGrid_pr = [min(kGrid_pdf),kGrid_pr(1)-1e6*eps,kGrid_pr];
        cdf = [0,0,cdf];
    elseif kGrid_pr(1)>min(kGrid_pdf)
        kGrid_pr = [min(kGrid_pdf),kGrid_pr];
        cdf = [0,cdf];
    end
    if kGrid_pr(end)<max(kGrid_pdf) && kGrid_pr(end)+1e6*eps<max(kGrid_pdf)
        kGrid_pr = [kGrid_pr,kGrid_pr(end)+1e6*eps,max(kGrid_pdf)];
        cdf = [cdf,1,1];
    elseif kGrid_pr(end)<max(kGrid_pdf)
        kGrid_pr = [kGrid_pr,max(kGrid_pdf)];
        cdf = [cdf,1];
    end
    cdf_cond_new(j,:) = interp1(kGrid_pr,cdf,kGrid_pdf);
    cdf_cond_new(j,end) = 1;
    pdf_cond_new(j,:) = [cdf_cond_new(j,1),cdf_cond_new(j,2:end)-cdf_cond_new(j,1:end-1)];       
end
end