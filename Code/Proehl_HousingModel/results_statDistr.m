% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This file computes the expected ergodic distributions and Euler equation 
% errors of te different algorithm solutions.
%__________________________________________________________________________
function results_statDistr(StaticParams,folder,order)

% truncation order 
M = order;

%--------------------------------------------------------------------------
% Computing the stationary distribution for our PPA and PFI with polynomial
% chaos
%--------------------------------------------------------------------------
for order_ct=0:M
% Initialization___________________________________________________________
Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'));
weights =  zeros(2,Poly.phi_N+1);
weights(:,1) = 2.5;
if order_ct>0
    weights(:,2) = 1;
end
Sol.cdf_cond = [];
pdf = reshape(Poly.total_pdf,[1,numel(Poly.total_pdf)]);
for j=1:2
    [Grid,idx_u] = sort(reshape(squeeze(sum(repmat(weights(j,:)',[1,size(Poly.total,2),size(Poly.total,3),size(Poly.total,4)]).*Poly.total(1:order_ct+1,:,:,:),1)),[1,numel(Poly.total_pdf)]));
    cdf = cumsum(pdf(idx_u));
    [Grid,idx_u] = unique(Grid,'last','legacy');
    cdf = cdf(idx_u);
    if Grid(1)>StaticParams.a_min
       Grid = [StaticParams.a_min,(Grid(1)-1e-10),Grid];
       cdf = [0,0,cdf];
    end
    if Grid(end)<StaticParams.a_max
       Grid = [Grid,(Grid(end)+1e-10),StaticParams.a_max];
       cdf = [cdf,1,1];
    end
    cdf = interp1(Grid,cdf,StaticParams.aGrid);
    Sol.cdf_cond = cummax(min(1,[Sol.cdf_cond;cdf;cdf]),2);
    Sol.cdf_cond(:,end)=1;
end
Sol.pdf_cond = [Sol.cdf_cond(:,1),Sol.cdf_cond(:,2:end)-Sol.cdf_cond(:,1:end-1)];

diffs_p = 1e8*ones(1,2);
iter = 1;
w_bounds = repmat([1e8;-1e8],[1,Poly.phi_N+1]);
steps = 1;

% Iteration________________________________________________________________
tstart = cputime;
while ((diffs_p(end,1) > StaticParams.criter_pdf)||(diffs_p(end,2)>1e-1)) && iter<2001
    a_prime = max(StaticParams.a_min,[getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,1,:)),Sol.idx,weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,2,:)),Sol.idx,weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,3,:)),Sol.idx,weights(2,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,4,:)),Sol.idx,weights(2,:))]);
    h_prime = max(StaticParams.a_min,[getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,1,:)),Sol.idx,weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2,:)),Sol.idx,weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,3,:)),Sol.idx,weights(2,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,4,:)),Sol.idx,weights(2,:))]);
        
    [~,~,pdf_cond_new,pdf_multistep] = calcStationaryDistr(a_prime,h_prime,StaticParams.P,...
                                     StaticParams.p,StaticParams,Sol.pdf_cond,steps);
    A_ss = [[StaticParams.ur(1),StaticParams.er(1)]*pdf_cond_new(1:2,:);...
            [StaticParams.ur(2),StaticParams.er(2)]*pdf_cond_new(3:4,:)]*StaticParams.aGrid';
    A = [[StaticParams.ur(1),StaticParams.er(1)]*pdf_multistep(1:2,:);...
         [StaticParams.ur(2),StaticParams.er(2)]*pdf_multistep(3:4,:)]*StaticParams.aGrid';
    diffs_p(iter,:) = [norm(pdf_multistep-Sol.pdf_cond,2),max(abs(A-A_ss))];
    if ((diffs_p(end,1) > StaticParams.criter_pdf)||(diffs_p(end,2)>1e-1))
        Sol.pdf_cond = 0.5.*Sol.pdf_cond+0.5.*pdf_multistep;
    else
        Sol.pdf_cond = pdf_cond_new;
    end

    cdf_new = cummax(min(1,cumsum(Sol.pdf_cond,2)),2);
    cdf_new(:,end) = 1;
    for i=1:2
        weights(i,:) = calcWeights(cdf_new(2*(i-1)+(1:2),:),Poly,StaticParams,i-1);
    end    
    w_bounds = [min(w_bounds(1,:),min(weights));...
                max(w_bounds(2,:),max(weights))];
    Sol.cdf_cond = cdf_new;

    iter = iter+1;
    if iter==2001
       Sol.DistrNoConvergence = {true,diffs_p(iter-1,:)}; 
    end
end

Sol.pdf = [[StaticParams.ur(1),StaticParams.er(1)]*Sol.pdf_cond(1:2,:);...
           [StaticParams.ur(2),StaticParams.er(2)]*Sol.pdf_cond(3:4,:)];
Sol.cdf_cond = cdf_new;
Sol.cdf = cummax(min(1,cumsum(Sol.pdf,2)),2);
Sol.cdf(1,Sol.cdf(1,:)==max(Sol.cdf(1,:))) = 1;
Sol.cdf(2,Sol.cdf(2,:)==max(Sol.cdf(2,:))) = 1;
Sol.Lorenz = Sol.cdf;
Sol.Lorenz(1,:)=(tril(repmat(Sol.pdf(1,:),[size(Sol.pdf,2),1]))*StaticParams.aGrid')'/(Sol.pdf(1,:)*StaticParams.aGrid');
Sol.Lorenz(2,:)=(tril(repmat(Sol.pdf(2,:),[size(Sol.pdf,2),1]))*StaticParams.aGrid')'/(Sol.pdf(2,:)*StaticParams.aGrid');
Sol.ergodic_gini = [sum(sum((Sol.pdf(1,:)'*Sol.pdf(1,:)).*abs(repmat(StaticParams.aGrid',[1,size(Sol.pdf,2)])-repmat(StaticParams.aGrid,[size(Sol.pdf,2),1])),1),2)...
            /(2*StaticParams.aGrid*Sol.pdf(1,:)');...
            sum(sum((Sol.pdf(2,:)'*Sol.pdf(2,:)).*abs(repmat(StaticParams.aGrid',[1,size(Sol.pdf,2)])-repmat(StaticParams.aGrid,[size(Sol.pdf,2),1])),1),2)...
            /(2*StaticParams.aGrid*Sol.pdf(2,:)')];

Sol.ergodic_weights = weights;
Sol.ergodic_r = [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,1)),Sol.idx,Sol.ergodic_weights(1,:)),...
             getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,2)),Sol.idx,Sol.ergodic_weights(2,:))];
Sol.ergodic_q_pr = [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.q_pr(:,1)),Sol.idx,Sol.ergodic_weights(1,:)),...
               getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.q_pr(:,2)),Sol.idx,Sol.ergodic_weights(2,:))];

c_approx = [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,Sol.ergodic_weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,Sol.ergodic_weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,Sol.ergodic_weights(2,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,Sol.ergodic_weights(2,:))];
c_approx = interp1(StaticParams.aGrid_pol,c_approx',StaticParams.aGrid)';
Sol.pdf_cond_c = calcEoPDistr(c_approx,Sol.cdf_cond,StaticParams.cGrid);
h_approx = [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,1,:)),Sol.idx,Sol.ergodic_weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2,:)),Sol.idx,Sol.ergodic_weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,3,:)),Sol.idx,Sol.ergodic_weights(2,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,4,:)),Sol.idx,Sol.ergodic_weights(2,:))];
h_approx = interp1(StaticParams.aGrid_pol,h_approx',StaticParams.aGrid)';
Sol.pdf_cond_h = calcEoPDistr(h_approx,Sol.cdf_cond,StaticParams.hGrid);
b_approx = [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.b(:,1,:)),Sol.idx,Sol.ergodic_weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.b(:,2,:)),Sol.idx,Sol.ergodic_weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.b(:,3,:)),Sol.idx,Sol.ergodic_weights(2,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.b(:,4,:)),Sol.idx,Sol.ergodic_weights(2,:))];
b_approx = interp1(StaticParams.aGrid_pol,b_approx',StaticParams.aGrid)';
Sol.pdf_cond_b = calcEoPDistr(b_approx,Sol.cdf_cond,StaticParams.bGrid);

Sol.time_distr = cputime-tstart;
clear a_prime c_prime pdf_cond_new cdf_new weights diffs iter ;

save(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'),'-struct','Sol');
end
end