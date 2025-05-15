% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function computes the new-period weights in the simulation for the 
% Euler equation errors.
%__________________________________________________________________________
function w = calcWeights(cdf,Poly,StaticParams,z_ag,kGrid)
if nargin>4
    w=zeros(2,Poly.phi_N+1);  
    pdf = calcEoPDistr(kGrid,cdf,StaticParams.kGrid);
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
        if i==1 || (i==2 && z_ag==0)
            cdf_aux(i,:,:) = cumsum(cdf_aux(i,:,:),3); 
            cdf_aux(i,2:end,:) = cdf_aux(i,2:end,:) + repmat(cumsum(cdf_aux(i,1:end-1,end)),[1,1,size(cdf_aux,3)]); 
        else
            cdf_aux(i,:,:) = cumsum(cdf_aux(i,:,:),2); 
            cdf_aux(i,:,2:end) = cdf_aux(i,:,2:end) + repmat(cumsum(cdf_aux(i,end,1:end-1)),[1,size(cdf_aux,3),1]); 
        end
    else
        cdf_aux(i,:) = cumsum(cdf_aux(i,:)); 
    end
    cdf_aux(i,:,:) = min(1,cdf_aux(i,:,:)./max(max(cdf_aux(i,:,:))));
end

Grid_new = zeros(size(Poly.total_pdf));
Grid_new(1,:,:,:) = projectDistr(StaticParams.kGrid,cdf(1,:),cdf_aux(1,:,:));
Grid_new(2,:,:,:) = projectDistr(StaticParams.kGrid,cdf(z_ag+1,:),cdf_aux(2,:,:));
Grid_new(3,:,:,:) = projectDistr(StaticParams.kGrid,cdf(2,:),cdf_aux(3,:,:));

w(1) = StaticParams.kGrid*[cdf(:,1),cdf(:,2:end)-cdf(:,1:end-1)]'*StaticParams.p(2*z_ag+(1:2))*2;
for idx=2:Poly.phi_N+1
    w(idx) = sum(sum(sum(Poly.total_pdf.*Grid_new.*squeeze(Poly.total(idx,:,:,:)),3),2),1)/sum(sum(sum(Poly.total_pdf.*squeeze(Poly.total(idx,:,:,:)).^2,3),2),1);
end
end