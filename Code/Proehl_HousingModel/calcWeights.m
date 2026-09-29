% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This function computes the new-period weights in the iteration for the 
% expected ergodic distribution and the simulation for the Euler equation
% errors.
%__________________________________________________________________________
function w = calcWeights(cdf,Poly,StaticParams,z_ag,aGrid)
if nargin>4
    w=zeros(2,Poly.phi_N+1);  
    pdf = calcEoPDistr(aGrid,cdf,StaticParams.aGrid);
    for j=1:2
        pdf_new = (repmat(1./(StaticParams.p(2*z_ag+(1:2))'*StaticParams.P(2*z_ag+(1:2),2*(j-1)+(1:2)))',[1,2])...
                       .*StaticParams.P(2*z_ag+(1:2),2*(j-1)+(1:2))')...
                       *(repmat(StaticParams.p(2*z_ag+(1:2)),[1,size(pdf,2)])...
                       .*squeeze(pdf));
        pdf_new = pdf_new.*(pdf_new>eps);
        pdf_new = pdf_new./repmat(sum(pdf_new,2),[1,size(pdf_new,2)]);
        cdf = cummax(min(1,cumsum(pdf_new,2)),2);
        cdf(1,cdf(1,:)==max(cdf(1,:))) = 1;
        cdf(2,cdf(2,:)==max(cdf(2,:))) = 1; 
        w(j,:) = calcOnce(cdf,j-1,Poly,StaticParams);
    end
else
    w = calcOnce(cdf,z_ag,Poly,StaticParams);
end
end

function w=calcOnce(cdf,z_ag,Poly,StaticParams)
w=zeros(1,Poly.phi_N+1);
cdf_aux = Poly.total_pdf;
for i=1:size(cdf_aux,1)
    if length(size(Poly.total_pdf))==3
        n_diags = 2*size(Poly.total_pdf,2)-1;
        idx = [];
        for j=1:n_diags
            a = max(1,size(Poly.total_pdf,2)-n_diags+j);
            b = min(size(Poly.total_pdf,2),j);
            idx = [idx;sub2ind([size(Poly.total_pdf,2),size(Poly.total_pdf,2)],a:b,j+1-(a:b))'];
        end
        aux = squeeze(cdf_aux(i,:,:));
        c_aux = cumsum(aux(idx));
        aux(idx) = c_aux;
        cdf_aux(i,:,:) = aux;
    else
        cdf_aux(i,:) = cumsum(cdf_aux(i,:)); 
    end
    cdf_aux(i,:,:) = min(1,cdf_aux(i,:,:)./max(max(cdf_aux(i,:,:))));
end

Grid_new = zeros(size(Poly.total_pdf));
Grid_new(1,:,:,:) = projectDistr(StaticParams.aGrid,cdf(1,:),cdf_aux(1,:,:));
Grid_new(2,:,:,:) = projectDistr(StaticParams.aGrid,cdf(z_ag+1,:),cdf_aux(2,:,:));
Grid_new(3,:,:,:) = projectDistr(StaticParams.aGrid,cdf(2,:),cdf_aux(3,:,:));

w(1) = StaticParams.aGrid*[cdf(:,1),cdf(:,2:end)-cdf(:,1:end-1)]'*StaticParams.p(2*z_ag+(1:2))*2;
for idx=2:Poly.phi_N+1
    w(idx) = sum(sum(sum(Poly.total_pdf.*Grid_new.*squeeze(Poly.total(idx,:,:,:)),3),2),1)/sum(sum(sum(Poly.total_pdf.*squeeze(Poly.total(idx,:,:,:)).^2,3),2),1);
end
end