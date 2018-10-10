% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This file implements the policy function iteration for the model with
% aggregate shocks.
%__________________________________________________________________________
function [c_new,k_prime_new,p_new,diffs] ...
         = OptPol_AggShock_PFI(pr,c,k_prime,basic_grid,StaticParams,Poly,distrGrid,idx,wealth,wage)

diffs = 1e8.*ones(2,1);
iter=0;
p_new = pr;
c_pr = zeros([size(k_prime,1),4,4,size(k_prime,3)]);
k_prime_new = zeros(size(k_prime));
weights_new = zeros(size(k_prime,1),2,2,Poly.phi_N+1);

grid = linspace(eps,100,10001);    
u_der_grid = 1./(grid.^StaticParams.gamma);   
u_der = @(x) interp1(grid,u_der_grid,max(eps,x),'linear','extrap');
u_der_inv = @(x) interp1(u_der_grid,grid,max(eps,x),'linear','extrap');

options = optimoptions('fsolve','Display','off','FunctionTolerance',StaticParams.criter_pdf);
	
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
       c_new = squeeze(zeros(size(k_prime(1,:,:,:))));
       c_new(1,:,:) = StaticParams.beta...
                    .*sum(repmat(permute(StaticParams.P(1,:),[3,1,2]),[1,1,1,length(StaticParams.kGrid_pol)])...
                    .*u_der(c_pr(m,1,:,:,:)),3);
       c_new(2,:,:) = StaticParams.beta...
                    .*sum(repmat(permute(StaticParams.P(2,:),[3,1,2]),[1,1,1,length(StaticParams.kGrid_pol)])...
                    .*u_der(c_pr(m,2,:,:,:)),3);
       c_new(3,:,:) = StaticParams.beta...
                    .*sum(repmat(permute(StaticParams.P(3,:),[3,1,2]),[1,1,1,length(StaticParams.kGrid_pol)])...
                    .*u_der(c_pr(m,3,:,:,:)),3);
       c_new(4,:,:) = StaticParams.beta...
                    .*sum(repmat(permute(StaticParams.P(4,:),[3,1,2]),[1,1,1,length(StaticParams.kGrid_pol)])...
                    .*u_der(c_pr(m,4,:,:,:)),3);
       c_prime_large = interp1(StaticParams.kGrid_pol,max(eps,u_der_inv(c_new))',StaticParams.kGrid)';   
       C = calcAggregates(c_prime_large,squeeze(basic_grid(m,:,:,:,:)),StaticParams,Poly);
       E = [2*StaticParams.p(1:2)'*squeeze(wage(m,1:2,1))';2*StaticParams.p(3:4)'*squeeze(wage(m,3:4,1))'];
       p_new(m,:) = fsolve(@(p) solveEquCond(squeeze(wealth(m,:,:)),c_new,p,squeeze(basic_grid(m,:,:,:,:)),StaticParams,Poly),...
                    (E./C).^StaticParams.gamma,options);
       [~,k_prime_new(m,:,:)] = solveEquCond(squeeze(wealth(m,:,:)),c_new,p_new(m,:),squeeze(basic_grid(m,:,:,:,:)),StaticParams,Poly);
   end
   c_new = max(eps,wealth-k_prime_new.*repmat(p_new(:,[1 1 2 2]),[1,1,length(StaticParams.kGrid_pol)]));
   
   iter = iter+1;
   diffs(iter,:) = max(max(max(max(abs(c-c_new)))));
   
   c = 0.4*c +0.6*c_new;
   k_prime = 0.4*k_prime +0.6*k_prime_new;
end
end

function [A,k_prime_new] = solveEquCond(wealth,c_new,p,basic_grid,StaticParams,Poly)
   grid = linspace(eps,100,10001);    
   u_der_grid = 1./(grid.^StaticParams.gamma);   
   u_der_inv = @(x) interp1(u_der_grid,grid,max(eps,x),'linear','extrap');
    
   k_prime_new = [max(StaticParams.k_min,(wealth(1:2,:)-max(eps,u_der_inv(c_new(1:2,:,:)./p(1))))./p(1));...
                  max(StaticParams.k_min,(wealth(3:4,:)-max(eps,u_der_inv(c_new(3:4,:,:)./p(2))))./p(2))];
   k_prime_large = interp1(StaticParams.kGrid_pol,k_prime_new',StaticParams.kGrid)';   
   A = calcAggregates(k_prime_large,basic_grid,StaticParams,Poly);
end