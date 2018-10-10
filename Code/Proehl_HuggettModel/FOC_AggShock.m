% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pr�hl
%
% AUTHOR Elisabeth Pr�hl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This function computes the first-order conditions to solve for in the
% proximal point algorithm with aggregate shocks.
%__________________________________________________________________________
function [f,J] = FOC_AggShock(sol,k_prime,c,p,y,y2,agK_m,grid,m,StaticParams,Sol,Poly,c_pr_der)

sol_k = reshape(sol,size(k_prime'))';
sol_k_large = interp1(StaticParams.kGrid_pol,sol_k',StaticParams.kGrid)';   
[~,weights_new] = calcAggregates(sol_k_large,grid,StaticParams,Poly);
weights_new(:,:,1) = 0;

c_prime = @(z_idx,z,k,w)...
           interp1(StaticParams.kGrid_pol,...
           getCApprox(Sol.distrGrid,squeeze(c(:,2*z(1)+z(2)+1,:)),...
           Sol.idx,squeeze(w((z_idx>2)+1,z(1)+1,:))'),...
           max(StaticParams.k_min,k),'linear','extrap');
c_pr = zeros([4,size(c,2),size(c,3)]);
for z_idx=1:4
    c_pr(z_idx,:,:) = [c_prime(z_idx,[0,0],sol_k(z_idx,:),weights_new);...
                       c_prime(z_idx,[0,1],sol_k(z_idx,:),weights_new);...
                       c_prime(z_idx,[1,0],sol_k(z_idx,:),weights_new);...
                       c_prime(z_idx,[1,1],sol_k(z_idx,:),weights_new)] ;
end   

f_foc_k = @(z_idx)...
         +p((z_idx>2)+1)*max(eps,squeeze(Sol.wealth(agK_m,z_idx,:))'-p((z_idx>2)+1)*sol_k(z_idx,:)).^(-StaticParams.gamma)...
         -StaticParams.beta*(StaticParams.P(z_idx,:)*...
          (squeeze(c_pr(z_idx,:,:)).^(-StaticParams.gamma)))...
         +1/m*(sol_k(z_idx,:)-k_prime(z_idx,:))...
         +((sol_k(z_idx,:)-StaticParams.k_min)<=y(z_idx,:)/m).*(-y(z_idx,:)+m*(sol_k(z_idx,:)-StaticParams.k_min))...
         +((squeeze(Sol.wealth(agK_m,z_idx,:))'-p((z_idx>2)+1)*sol_k(z_idx,:))<=y2(z_idx,:)/m)...
          .*(-y2(z_idx,:)+m*(squeeze(Sol.wealth(agK_m,z_idx,:))'-p((z_idx>2)+1)*sol_k(z_idx,:)));

% function value of the Euler equation       
foc_k = [f_foc_k(1);f_foc_k(2);f_foc_k(3);f_foc_k(4)].*squeeze(abs(Sol.wealth(agK_m,:,:)));
f = reshape(foc_k',1,numel(foc_k))';

if nargout>1    
    f_foc_k_dk = @(z_idx)...
        +(StaticParams.gamma)*p((z_idx>2)+1)^2.*(squeeze(Sol.wealth(agK_m,z_idx,:))'-p((z_idx>2)+1)*sol_k(z_idx,:)>eps)...
         .*max(eps,squeeze(Sol.wealth(agK_m,z_idx,:))'-p((z_idx>2)+1)*sol_k(z_idx,:)).^(-StaticParams.gamma-1)...
        -StaticParams.beta*(StaticParams.P(z_idx,:)*...
         ((-StaticParams.gamma).*squeeze(c_pr(z_idx,:,:)).^(-StaticParams.gamma-1)...
         .*squeeze(c_pr_der(z_idx,:,:))))...
        +1/m+((sol_k(z_idx,:)-StaticParams.k_min)<=y(z_idx,:)/m).*m...
        +((squeeze(Sol.wealth(agK_m,z_idx,:))'-p((z_idx>2)+1)*sol_k(z_idx,:))<=y2(z_idx,:)/m).*(-p((z_idx>2)+1)*m);...
      
    % Jacobian of the Euler equation
    foc_k_dk = [f_foc_k_dk(1);f_foc_k_dk(2);f_foc_k_dk(3);f_foc_k_dk(4)].*squeeze(abs(Sol.wealth(agK_m,:,:)));
    J = sparse(1:length(sol),1:length(sol),reshape(foc_k_dk',numel(foc_k_dk),1));
end
end