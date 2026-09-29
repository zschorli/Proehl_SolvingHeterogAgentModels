% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This function calculates the law of motion from the beginning-of-period
% distribution to the end-of-period distribution after the agents apply 
% their policy.
%__________________________________________________________________________

function pdf_cond_new = calcEoPDistr(a_prime,cdf_cond_BoP,aGrid_pdf)

cdf_cond_new = zeros(size(aGrid_pdf));
pdf_cond_new = zeros(size(aGrid_pdf));

pdf_cond_aux = [cdf_cond_BoP(:,1),cdf_cond_BoP(:,2:end)-cdf_cond_BoP(:,1:end-1)]; 
[a_prime,idx] = sort(a_prime,2,"ascend");
idx = sub2ind(size(idx),repmat((1:size(idx,1))',[1,size(idx,2)]),idx);
pdf_cond_aux = pdf_cond_aux(idx);
cdf_cond_aux = cumsum(pdf_cond_aux,2);
cdf_cond_aux(:,end) = 1;

for j=1:size(a_prime,1)
    [aGrid_pr,idx] = unique(a_prime(j,:),'last','legacy');
    cdf = cdf_cond_aux(j,idx);
    if aGrid_pr(1)>min(aGrid_pdf) && aGrid_pr(1)-1e6*eps>min(aGrid_pdf) 
        aGrid_pr = [min(aGrid_pdf),aGrid_pr(1)-1e6*eps,aGrid_pr];
        cdf = [0,0,cdf];
    elseif aGrid_pr(1)>min(aGrid_pdf)
        aGrid_pr = [min(aGrid_pdf),aGrid_pr];
        cdf = [0,cdf];
    end
    if aGrid_pr(end)<max(aGrid_pdf) && aGrid_pr(end)+1e6*eps<max(aGrid_pdf)
        aGrid_pr = [aGrid_pr,aGrid_pr(end)+1e6*eps,max(aGrid_pdf)];
        cdf = [cdf,1,1];
    elseif aGrid_pr(end)<max(aGrid_pdf)
        aGrid_pr = [aGrid_pr,max(aGrid_pdf)];
        cdf = [cdf,1];
    end
    cdf_cond_new(j,:) = interp1(aGrid_pr,cdf,aGrid_pdf);
    cdf_cond_new(j,end) = 1;
    pdf_cond_new(j,:) = [cdf_cond_new(j,1),cdf_cond_new(j,2:end)-cdf_cond_new(j,1:end-1)];       
end
end