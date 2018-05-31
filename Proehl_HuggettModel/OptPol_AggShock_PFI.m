% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This file implements the policy function iteration for the model with
% aggregate shocks.
%__________________________________________________________________________
function [c_new,k_prime_new,p_new,diffs] ...
         = OptPol_AggShock_PFI(pr,c,k_prime,basic_grid,StaticParams,Poly,distrGrid,idx,wealth)

diffs = 1e8.*ones(2,1);
iter=0;
p_new = pr;
c_pr = zeros([size(k_prime,1),4,4,size(k_prime,3)]);
k_prime_new = zeros(size(k_prime));
weights_new = zeros(size(k_prime,1),2,2,Poly.phi_N+1);
	
while max(diffs(max(iter,1),:)) > StaticParams.criter_k
   parfor m=1:size(k_prime,1)
       k_prime_large = interp1(StaticParams.kGrid_pol,squeeze(k_prime(m,:,:))',StaticParams.kGrid)';   
       [~,w] = calcAggregates(k_prime_large,squeeze(basic_grid(m,:,:,:,:)),StaticParams,Poly);
       w(:,:,1) = 0;
       weights_new(m,:,:,:) = w;
       c_prime = @(z_idx,z,k,w)...
           interp1(StaticParams.kGrid_pol,...
           getCApprox(distrGrid,squeeze(c(:,2*z(1)+z(2)+1,:)),...
           idx,squeeze(w((z_idx>2)+1,z(1)+1,:))'),...
           max(StaticParams.k_min,k),'linear','extrap');
       for z_idx=1:4
           c_pr(m,z_idx,:,:) = max(eps,[c_prime(z_idx,[0,0],squeeze(k_prime(m,z_idx,:))',squeeze(weights_new(m,:,:,:)));...
                                  c_prime(z_idx,[0,1],squeeze(k_prime(m,z_idx,:))',squeeze(weights_new(m,:,:,:)));...
                                  c_prime(z_idx,[1,0],squeeze(k_prime(m,z_idx,:))',squeeze(weights_new(m,:,:,:)));...
                                  c_prime(z_idx,[1,1],squeeze(k_prime(m,z_idx,:))',squeeze(weights_new(m,:,:,:)))]);
       end
   end
   parfor m=1:size(k_prime,1)
       % outer loop to find the price
       options = optimset('Display','off','Jacobian','off','TolFun',eps,'TolX',eps,'MaxIter',400);
       price = p_new(m,:);
       pr_new1=fzero(@(p) solveFOC_PFI(p,m,0,k_prime,c_pr,StaticParams,Poly,basic_grid,wealth),...
                      price(1),options);
       pr_new2=fzero(@(p) solveFOC_PFI(p,m,1,k_prime,c_pr,StaticParams,Poly,basic_grid,wealth),...
                      price(2),options); 
       p_new(m,:) = [pr_new1,pr_new2];
       [~,sol1] = solveFOC_PFI(pr_new1,m,0,k_prime,c_pr,StaticParams,Poly,basic_grid,wealth);
       [~,sol2] = solveFOC_PFI(pr_new2,m,1,k_prime,c_pr,StaticParams,Poly,basic_grid,wealth);
       k_prime_new(m,:,:) = [sol1;sol2];
   end
   c_new = max(eps,wealth-k_prime_new.*repmat(p_new(:,[1 1 2 2]),[1,1,length(StaticParams.kGrid_pol)]));
   
   iter = iter+1;
   diffs(iter,:) = max(max(max(max(abs(c-c_new)))));
   
   c = 0.4*c +0.6*c_new;
   k_prime = 0.4*k_prime +0.6*k_prime_new;
end
end

function [A,sol]=solveFOC_PFI(p,m,z_ag,k_prime,c_pr,StaticParams,Poly,basic_grid,wealth)
% inner loop to find the policy
grid = linspace(eps,100,10001);    
u_der_grid = 1./(grid.^StaticParams.gamma);   
u_der = @(x) interp1(grid,u_der_grid,max(eps,x),'linear','extrap');
u_der_inv = @(x) interp1(u_der_grid,grid,max(eps,x),'linear','extrap');

c_new = squeeze(zeros(size(k_prime(1,1:2,:,:))));
c_new(1,:,:) = max(eps,u_der_inv(StaticParams.beta./p...
            .*sum(repmat(permute(StaticParams.P(z_ag*2+1,:),[3,1,2]),[1,1,1,length(StaticParams.kGrid_pol)])...
            .*u_der(c_pr(m,z_ag*2+1,:,:,:)),3)));
c_new(2,:,:) = max(eps,u_der_inv(StaticParams.beta./p...
            .*sum(repmat(permute(StaticParams.P(z_ag*2+2,:),[3,1,2]),[1,1,1,length(StaticParams.kGrid_pol)])...
            .*u_der(c_pr(m,z_ag*2+2,:,:,:)),3)));
sol = max(StaticParams.k_min,(squeeze(wealth(m,z_ag*2+(1:2),:))-c_new)./p);
sol_k_large = interp1(StaticParams.kGrid_pol,[sol;sol]',StaticParams.kGrid)';   
K = calcAggregates(sol_k_large,squeeze(basic_grid(m,:,:,:)),StaticParams,Poly);        
A=K(z_ag+1);
end