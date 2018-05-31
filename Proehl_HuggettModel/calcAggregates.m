% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This function computes the next-period aggregate capital and weights.
%__________________________________________________________________________
function [K,w,K1,K2] = calcAggregates(sol_k,grid,StaticParams,Poly)    
K = zeros(2,size(sol_k,3)); % aggregate capital
K1 = zeros(2,size(sol_k,3)); % aggregate capital of unemployed
K2 = zeros(2,size(sol_k,3)); % aggregate capital of employed
if nargout>1
    w = zeros(2,2,Poly.phi_N+1,size(sol_k,3));
end
if size(sol_k,1)<=2
   sol_k = repmat(sol_k,[2,1,1]); 
end
pdf = [reshape(Poly.total_pdf(1,:,:,:),[1,numel(Poly.total_pdf(1,:,:,:))]);...
       reshape(Poly.total_pdf(2,:,:,:),[1,numel(Poly.total_pdf(1,:,:,:))]);...
       reshape(Poly.total_pdf(3,:,:,:),[1,numel(Poly.total_pdf(1,:,:,:))])];

for s=1:size(sol_k,3)
    k_prime = @(z_idx,grid) interp1(StaticParams.kGrid,squeeze(sol_k(z_idx,:,s))',max(StaticParams.k_min,grid),'linear','extrap');
    for j=1:2
        grid_new = [k_prime(2*(j-1)+1,reshape(grid(1,:,:,:),[1,numel(grid(1,:,:,:))]));...
                    k_prime(2*(j-1)+j,reshape(grid(2,:,:,:),[1,numel(grid(1,:,:,:))]));...
                    k_prime(2*(j-1)+2,reshape(grid(3,:,:,:),[1,numel(grid(1,:,:,:))]))];
        
        K(j,s) = sum(sum(grid_new.*pdf,2),1);
        K1(j,s) = sum(sum(grid_new(1:3-j,:).*pdf(1:3-j,:),2),1)/StaticParams.ur(j);
        K2(j,s) = sum(sum(grid_new(4-j:3,:).*pdf(4-j:3,:),2),1)/StaticParams.er(j);
        
        if nargout>1            
            for l=1:2
                grid_new = squeeze(Poly.T(j,l,:,:))'*(grid_new.*pdf);
                grid_new = grid_new./max(eps,squeeze(Poly.T(j,l,:,:))'*pdf);
                Grid_new = zeros(size(Poly.total_pdf));
                Grid_new(1,:,:,:) = reshape(grid_new(1,:),size(Poly.total_pdf(1,:,:,:)));
                Grid_new(2,:,:,:) = reshape(grid_new(2,:),size(Poly.total_pdf(1,:,:,:)));
                Grid_new(3,:,:,:) = reshape(grid_new(3,:),size(Poly.total_pdf(1,:,:,:)));
                
                for idx=1:Poly.phi_N+1
                    w(j,l,idx,s) = sum(sum(sum(sum(Poly.total_pdf.*Grid_new.*squeeze(Poly.total(idx,:,:,:,:)),4),3),2),1)/sum(sum(sum(sum(Poly.total_pdf.*squeeze(Poly.total(idx,:,:,:,:)).^2,4),3),2),1);
                end
            end
        end
    end
end
end