% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This file computes the expected ergodic distributions and Euler equation 
% errors of te different algorithm solutions.
%__________________________________________________________________________
function results_statDistr(StaticParams,folder)

% truncation order 
M = 4;

%--------------------------------------------------------------------------
% Computing the stationary distribution for our PPA and PFI with polynomial
% chaos
%--------------------------------------------------------------------------
% Initialization___________________________________________________________
for method = 1:-1:0
for order_ct=M%1:M
if method==0
    Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PPA.mat'));
    preSol = load(strcat(folder,'/preSol_PPA.mat'));
    Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PPA.mat'));
else
    Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
    preSol = load(strcat(folder,'/preSol_PFI.mat'));
    Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'));
end

weights =  zeros(2,Poly.phi_N+1);
for i=1:2
    weights(i,:) = calcWeights((preSol.cdf_cond_b+preSol.cdf_cond_g)/2,Poly,StaticParams,i-1);
end 
Sol.cdf_cond = [];
pdf = reshape(Poly.total_pdf,[1,numel(Poly.total_pdf)]);
for j=1:2
    [Grid,idx_u] = sort(reshape(repmat(Poly.pdf{1}',[1,size(Poly.total,3),size(Poly.total,4)]).*squeeze(sum(repmat(weights(j,:)',[1,size(Poly.total,2),size(Poly.total,3),size(Poly.total,4)]).*Poly.total(1:order_ct+1,:,:,:),1)),[1,numel(Poly.total_pdf)]));
    cdf = cumsum(pdf(idx_u));
    [Grid,idx_u] = unique(Grid,'last','legacy');
    cdf = cdf(idx_u);
    if Grid(1)>StaticParams.k_min
       Grid = [StaticParams.k_min,(Grid(1)-1e-10),Grid];
       cdf = [0,0,cdf];
    end
    if Grid(end)<StaticParams.k_max
       Grid = [Grid,(Grid(end)+1e-10),StaticParams.k_max];
       cdf = [cdf,1,1];
    end
    cdf = interp1(Grid,cdf,StaticParams.kGrid);
    Sol.cdf_cond = cummax(min(1,[Sol.cdf_cond;cdf;cdf]),2);
end
Sol.cdf_cond(1,Sol.cdf_cond(1,:)==max(Sol.cdf_cond(1,:))) = 1;
Sol.cdf_cond(2,Sol.cdf_cond(2,:)==max(Sol.cdf_cond(2,:))) = 1;
Sol.cdf_cond(3,Sol.cdf_cond(3,:)==max(Sol.cdf_cond(3,:))) = 1;
Sol.cdf_cond(4,Sol.cdf_cond(4,:)==max(Sol.cdf_cond(4,:))) = 1;
Sol.pdf_cond = [Sol.cdf_cond(:,1),Sol.cdf_cond(:,2:end)-Sol.cdf_cond(:,1:end-1)];

diffs_p = 1e8*ones(2,1);
iter = 1;
w_bounds = repmat([1e8;-1e8],[1,Poly.phi_N+1]);


% Iteration________________________________________________________________
tic;
while max(diffs_p(end)) > StaticParams.criter_pdf && iter<1001
    k_prime = max(StaticParams.k_min,[getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,1,:)),Sol.idx,[0,weights(1,2:end)]);...
                     getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2,:)),Sol.idx,[0,weights(1,2:end)]);...
                     getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,3,:)),Sol.idx,[0,weights(2,2:end)]);...
                     getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,4,:)),Sol.idx,[0,weights(2,2:end)])]);
    pdf_cond_new = calcEoPDistr(interp1(StaticParams.kGrid_pol,k_prime',StaticParams.kGrid)',...
                   Sol.cdf_cond,StaticParams.kGrid);
    pdf_cond_new = (repmat(1./(StaticParams.p'*StaticParams.P)',[1,4])...
               .*StaticParams.P')...
               *(repmat(StaticParams.p,[1,size(Sol.pdf_cond,2)])...
               .*pdf_cond_new);
    pdf_cond_new = pdf_cond_new.*(pdf_cond_new>eps);
    pdf_cond_new = pdf_cond_new./repmat(sum(pdf_cond_new,2),[1,size(pdf_cond_new,2)]);
    cdf_new = cummax(min(1,cumsum(pdf_cond_new,2)),2);
    cdf_new(1,cdf_new(1,:)==max(cdf_new(1,:))) = 1;
    cdf_new(2,cdf_new(2,:)==max(cdf_new(2,:))) = 1;
    cdf_new(3,cdf_new(3,:)==max(cdf_new(3,:))) = 1;
    cdf_new(4,cdf_new(4,:)==max(cdf_new(4,:))) = 1;
    for i=1:2
        weights(i,:) = calcWeights(cdf_new(2*(i-1)+(1:2),:),Poly,StaticParams,i-1);
    end    
    w_bounds = [min(w_bounds(1,:),min(weights));...
                max(w_bounds(2,:),max(weights))];

    diffs_p(iter) = norm(pdf_cond_new-Sol.pdf_cond,2);
    Sol.pdf_cond = pdf_cond_new;
    Sol.cdf_cond = cdf_new;
    iter = iter+1;
    if iter==1001
       Sol.DistrNoConvergence = {true,diffs_p}; 
    end
end
Sol.pdf = [[StaticParams.ur(1),StaticParams.er(1)]*Sol.pdf_cond(1:2,:);...
           [StaticParams.ur(2),StaticParams.er(2)]*Sol.pdf_cond(3:4,:)];
Sol.cdf_cond = cdf_new;
Sol.cdf = cummax(min(1,cumsum(Sol.pdf,2)),2);
Sol.cdf(1,Sol.cdf(1,:)==max(Sol.cdf(1,:))) = 1;
Sol.cdf(2,Sol.cdf(2,:)==max(Sol.cdf(2,:))) = 1;
Sol.statPrices = [getCApprox(Sol.distrGrid,squeeze(Sol.p(:,1)),Sol.idx,[0,weights(1,2:end)]),...
                  getCApprox(Sol.distrGrid,squeeze(Sol.p(:,2)),Sol.idx,[0,weights(2,2:end)])];
Sol.time_distr = toc;
clear k_prime c_prime pdf_cond_new cdf_new weights diffs iter ;

if method==0
    save(strcat(folder,'/Sol',num2str(order_ct),'_PPA.mat'),'-struct','Sol');
else
    save(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'),'-struct','Sol');
end
end
end

clear cdf1 cdf2 cdf_new1 cdf_new2 diffs_p unemplGrid unemplGrid_new ...
      emplGrid emplGrid_new i j idx_cdf1 idx_cdf2 idx_s1 idx_s2 k_e k_u ...
      max_weights min_weights order_ct cdf_e cdf_u idx_e idx_u method; 

end