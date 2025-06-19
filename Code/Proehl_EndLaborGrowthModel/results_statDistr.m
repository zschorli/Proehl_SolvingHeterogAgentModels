% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file computes the expected ergodic distributions of the different 
% algorithm solutions.
%__________________________________________________________________________
function results_statDistr(StaticParams,case_nr,folder,order)

% truncation order 
M = order;

%--------------------------------------------------------------------------
% Computing the stationary distribution for our PPA and PFI with polynomial
% chaos
%--------------------------------------------------------------------------
% Initialization___________________________________________________________
if case_nr==1
    start = 0;
else
    start = 1;
end
for order_ct=start:M
Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'));
if order_ct==0
    switch case_nr
    case {1,5,6,14,19,20,21,29,30,31}
        count = 1;
    case {2,3,4,7}
        count = case_nr;
    otherwise
        count = 7;
    end
    preSol =load(strcat('res_case',num2str(count),'/preSol_PFI.mat'));
    weights =  [StaticParams.kGrid*preSol.pdf_cond_b'*[StaticParams.ur(1);StaticParams.er(1)];...
                StaticParams.kGrid*preSol.pdf_cond_g'*[StaticParams.ur(2);StaticParams.er(2)]];
else
    if order_ct>1
        Solprev = load(strcat(folder,'/Sol',num2str(order_ct-1),'_PFI.mat'));
        weights = [calcWeights(Solprev.cdf_cond(1:2,:),Poly,StaticParams,0);...
                   calcWeights(Solprev.cdf_cond(3:4,:),Poly,StaticParams,1)];
    else 
        switch case_nr
        case {1,5,6,14,19,20,21,29,30,31}
            count = 1;
        case {2,3,4,7}
            count = case_nr;
        otherwise
            count = 7;
        end
        preSol =load(strcat('res_case',num2str(count),'/preSol_PFI.mat'));
        weights =  [StaticParams.kGrid*preSol.pdf_cond_b'*[StaticParams.ur(1);StaticParams.er(1)],1,zeros(1,order_ct-1);...
                    StaticParams.kGrid*preSol.pdf_cond_g'*[StaticParams.ur(2);StaticParams.er(2)],1,zeros(1,order_ct-1)];
    end
end
Sol.cdf_cond = [];
for j=1:2
    [Grid,idx_u] = sort(Poly.pdf{1}*squeeze(sum(repmat(weights(j,:)',[1,size(Poly.total,2),size(Poly.total,3)]).*Poly.total(1:order_ct+1,:,:),1)));
    cdf = cumsum(Poly.pdf{2}(idx_u));
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

diffs_p = 1e8*ones(2,2);
iter = 1;
w_bounds = repmat([1e8;-1e8],[1,Poly.phi_N+1]);
steps = 10;

% Iteration________________________________________________________________
tic;
while (diffs_p(end,1) > StaticParams.criter_pdf) && (iter<151)
    k_prime = [max(StaticParams.k_min,getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,1,:)),Sol.idx,weights(1,:)));...
               max(StaticParams.k_min,getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2,:)),Sol.idx,weights(1,:)));...
               max(StaticParams.k_min2,getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,3,:)),Sol.idx,weights(2,:)));...
               max(StaticParams.k_min2,getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,4,:)),Sol.idx,weights(2,:)))];
    l_aux = max(0,[getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,1,:)),Sol.idx,weights(1,:));...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,2,:)),Sol.idx,weights(1,:));...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,3,:)),Sol.idx,weights(2,:));...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,4,:)),Sol.idx,weights(2,:))]);
    [~,~,pdf_cond_new,pdf_multistep] = calcStationaryDistr(k_prime,l_aux,StaticParams.P,...
                                     StaticParams.p,StaticParams,Sol.pdf_cond,steps);

    K_ss = [[StaticParams.ur(1),StaticParams.er(1)]*pdf_cond_new(1:2,:);...
            [StaticParams.ur(2),StaticParams.er(2)]*pdf_cond_new(3:4,:)]*StaticParams.kGrid';
    K = [[StaticParams.ur(1),StaticParams.er(1)]*pdf_multistep(1:2,:);...
         [StaticParams.ur(2),StaticParams.er(2)]*pdf_multistep(3:4,:)]*StaticParams.kGrid';
    diffs_p(iter,:) = [norm(pdf_multistep-Sol.pdf_cond,2),max(abs(K-K_ss))];
    if ((diffs_p(end,1) > StaticParams.criter_pdf)||(diffs_p(end,2)>1e-1))
        Sol.pdf_cond = pdf_multistep;
    else
        Sol.pdf_cond = pdf_cond_new;
    end

    cdf_new = cummax(min(1,cumsum(Sol.pdf_cond,2)),2);
    cdf_new(1,cdf_new(1,:)==max(cdf_new(1,:))) = 1;
    cdf_new(2,cdf_new(2,:)==max(cdf_new(2,:))) = 1;
    cdf_new(3,cdf_new(3,:)==max(cdf_new(3,:))) = 1;
    cdf_new(4,cdf_new(4,:)==max(cdf_new(4,:))) = 1;
    for i=1:2
        weights(i,:) = calcWeights(cdf_new(2*(i-1)+(1:2),:),Poly,StaticParams,i-1);
    end    
    w_bounds = [min(w_bounds(1,:),min(weights));...
                max(w_bounds(2,:),max(weights))];
    Sol.cdf_cond = cdf_new;
    iter = iter+1;
end
if (iter==152) || (diffs_p(end,2)>=20)
   Sol.DistrNoConvergence = {true,diffs_p(iter-1,:)}; 
end
[Yz,Yk] = ndgrid(1:4,StaticParams.kGrid);
l_aux = interpn(1:4,StaticParams.kGrid_pol,l_aux,Yz,Yk);
Sol.pdf_cond_l = calcEoPDistr(1-l_aux,Sol.cdf_cond,StaticParams.lGrid);
Sol.pdf_cond_l = Sol.pdf_cond_l(:,size(Sol.pdf_cond_l,2):-1:1);
Sol.cdf_cond_l = cummax(min(1,cumsum(Sol.pdf_cond_l,2)),2);
Sol.cdf_cond_l(1,Sol.cdf_cond_l(1,:)==max(Sol.cdf_cond_l(1,:))) = 1;
Sol.cdf_cond_l(2,Sol.cdf_cond_l(2,:)==max(Sol.cdf_cond_l(2,:))) = 1;
Sol.cdf_cond_l(3,Sol.cdf_cond_l(3,:)==max(Sol.cdf_cond_l(3,:))) = 1;
Sol.cdf_cond_l(4,Sol.cdf_cond_l(4,:)==max(Sol.cdf_cond_l(4,:))) = 1;
Sol.pdf = [[StaticParams.ur(1),StaticParams.er(1)]*Sol.pdf_cond(1:2,:);...
           [StaticParams.ur(2),StaticParams.er(2)]*Sol.pdf_cond(3:4,:)];
Sol.cdf = cummax(min(1,cumsum(Sol.pdf,2)),2);
Sol.cdf(1,Sol.cdf(1,:)==max(Sol.cdf(1,:))) = 1;
Sol.cdf(2,Sol.cdf(2,:)==max(Sol.cdf(2,:))) = 1;

Sol.time_distr = toc;
display(w_bounds);
save(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'),'-struct','Sol');
end
end