% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file computes the expected ergodic distributions and Euler equation 
% errors of te different algorithm solutions.
%__________________________________________________________________________
function results_statDistr(case_nr,StaticParams,folder)

% truncation order 
% switch case_nr
%     case 1
%         M = 3;
%         M_KS = 4;
%     case 3
%         M = 3;
%         M_KS = 4;
%     otherwise
        M = 3;
        M_KS = 4;
% end

% --------------------------------------------------------------------------
% Computing the stationary distribution for the Krusell-Smith Algo from
% Maliar, Maliar, Valli (2010)
% --------------------------------------------------------------------------
for order=1:M_KS
% Initialization___________________________________________________________
folder_KS = strcat('C:\Users\Elisabeth Proehl\Documents\GitHub\Proehl_SolvingHeterogAgentModels\Code\KrusellSmithByMaliarMaliarValli_',num2str(order),'mom');
KS = load(strcat(folder_KS,'/KS',num2str(order),'_Sol',num2str(case_nr),'.mat'));
KS.pdf_cond = zeros(4,length(StaticParams.kGrid));  
KS.pdf_cond([1,2],getCapIdx(StaticParams.kss_b,StaticParams.kGrid)) = 1;
KS.pdf_cond([3,4],getCapIdx(StaticParams.kss_g,StaticParams.kGrid)) = 1;
KS.cdf_cond = cumsum(KS.pdf_cond,2);
KS.cdf_cond(:,end) = 1;
KS.pdf = [[StaticParams.ur(1),StaticParams.er(1)]*KS.pdf_cond(1:2,:);...
          [StaticParams.ur(2),StaticParams.er(2)]*KS.pdf_cond(3:4,:)];
K = KS.pdf*StaticParams.kGrid';
diffs_p = 1e8*ones(2,2);
iter = 1;   
steps = 50;

% Iteration________________________________________________________________
while ((diffs_p(end,1) > StaticParams.criter_pdf)||(diffs_p(end,2)>1e-1)) && iter<2001
    switch order
        case 1
        k_prime = max(0,[interpn(KS.k,KS.km,KS.kprime(:,:,1,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),'spline')';...
                     interpn(KS.k,KS.km,KS.kprime(:,:,1,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),'spline')';...
                     interpn(KS.k,KS.km,KS.kprime(:,:,2,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),'spline')';...
                     interpn(KS.k,KS.km,KS.kprime(:,:,2,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),'spline')']);
        case 2
        sd = [sqrt(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^2);sqrt(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^2)];
        k_prime = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,1,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),'spline')';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,1,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),'spline')';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,2,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),'spline')';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,2,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),'spline')']);
        case 3
        sd = [sqrt(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^2);sqrt(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^2)];
        skew = [(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^3)/sd(1)^3;(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^3)/sd(2)^3];
        k_prime = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,1,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),min(KS.kskew_max,max(KS.kskew_min,skew(1))),'spline')';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,1,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),min(KS.kskew_max,max(KS.kskew_min,skew(1))),'spline')';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,2,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),min(KS.kskew_max,max(KS.kskew_min,skew(2))),'spline')';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,2,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),min(KS.kskew_max,max(KS.kskew_min,skew(2))),'spline')']);
        case 4
        sd = [sqrt(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^2);sqrt(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^2)];
        skew = [(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^3)/sd(1)^3;(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^3)/sd(2)^3];
        kurt = [(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^4)/sd(1)^4;(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^4)/sd(2)^4];
        k_prime = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,1,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),min(KS.kskew_max,max(KS.kskew_min,skew(1))),min(KS.kkurt_max,max(KS.kkurt_min,kurt(1))),'spline')';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,1,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),min(KS.kskew_max,max(KS.kskew_min,skew(1))),min(KS.kkurt_max,max(KS.kkurt_min,kurt(1))),'spline')';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,2,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),min(KS.kskew_max,max(KS.kskew_min,skew(2))),min(KS.kkurt_max,max(KS.kkurt_min,kurt(2))),'spline')';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,2,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),min(KS.kskew_max,max(KS.kskew_min,skew(2))),min(KS.kkurt_max,max(KS.kkurt_min,kurt(2))),'spline')']);
    end           
    [~,pdf_cond_new,pdf_multistep] = calcStationaryDistr(k_prime,StaticParams.P,...
                                     StaticParams.p,StaticParams,KS.pdf_cond,steps);

    K_ss = [[StaticParams.ur(1),StaticParams.er(1)]*pdf_cond_new(1:2,:);...
            [StaticParams.ur(2),StaticParams.er(2)]*pdf_cond_new(3:4,:)]*StaticParams.kGrid';
    K = [[StaticParams.ur(1),StaticParams.er(1)]*pdf_multistep(1:2,:);...
         [StaticParams.ur(2),StaticParams.er(2)]*pdf_multistep(3:4,:)]*StaticParams.kGrid';
    diffs_p(iter,:) = [norm(pdf_multistep-KS.pdf_cond,2),max(abs(K-K_ss))];
    if ((diffs_p(end,1) > StaticParams.criter_pdf)||(diffs_p(end,2)>1e-1)) && iter<2001
        KS.pdf_cond = pdf_multistep;
    else
        KS.pdf_cond = pdf_cond_new;
    end

    KS.pdf = [[StaticParams.ur(1),StaticParams.er(1)]*KS.pdf_cond(1:2,:);...
              [StaticParams.ur(2),StaticParams.er(2)]*KS.pdf_cond(3:4,:)];
    KS.cdf_cond = cummax(min(1,cumsum(KS.pdf_cond,2)),2);
    KS.cdf_cond(1,KS.cdf_cond(1,:)==max(KS.cdf_cond(1,:))) = 1;
    KS.cdf_cond(2,KS.cdf_cond(2,:)==max(KS.cdf_cond(2,:))) = 1;
    iter = iter+1;
    if iter==2001
       KS.DistrNoConvergence = {true,diffs_p(iter-1,:)}; 
    end
end
KS.cdf = cumsum(KS.pdf,2);
KS.cdf(:,end) = 1;
save(strcat(folder,'/KS',num2str(order),'_Sol',num2str(case_nr),'.mat'),'-struct','KS');
end

%--------------------------------------------------------------------------
% Computing the stationary distribution for our PPA and PFI with polynomial
% chaos
%--------------------------------------------------------------------------
% Initialization___________________________________________________________
for method = 1:-1:0
for order_ct=0:M
if method==0
    Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PPA.mat'));
    preSol = load(strcat(folder,'/preSol_PPA.mat'));
    Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PPA.mat'));
else
    Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
    preSol = load(strcat(folder,'/preSol_PFI.mat'));
    Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'));
end

if order_ct>0
    weights =  [StaticParams.kGrid*preSol.pdf_cond_b'*[StaticParams.ur(1);StaticParams.er(1)],1,zeros(1,Poly.phi_N-1);...
                StaticParams.kGrid*preSol.pdf_cond_g'*[StaticParams.ur(2);StaticParams.er(2)],1,zeros(1,Poly.phi_N-1)];
    % weights =  [StaticParams.kGrid*preSol.pdf_cond_b'*[StaticParams.ur(1);StaticParams.er(1)],2.5,0.04,zeros(1,Poly.phi_N-2);...
    %             StaticParams.kGrid*preSol.pdf_cond_g'*[StaticParams.ur(2);StaticParams.er(2)],2.5,0.04,zeros(1,Poly.phi_N-2)];
else
    weights =  [StaticParams.kGrid*preSol.pdf_cond_b'*[StaticParams.ur(1);StaticParams.er(1)];...
                StaticParams.kGrid*preSol.pdf_cond_g'*[StaticParams.ur(2);StaticParams.er(2)]];
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
steps = 50;

% Iteration________________________________________________________________
tic;
while ((diffs_p(end,1) > StaticParams.criter_pdf)||(diffs_p(end,2)>1e-1)) && iter<2001
    k_prime = max(0,[getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,1,:)),Sol.idx,weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2,:)),Sol.idx,weights(1,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,3,:)),Sol.idx,weights(2,:));...
                     getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,4,:)),Sol.idx,weights(2,:))]);
    [~,pdf_cond_new,pdf_multistep] = calcStationaryDistr(k_prime,StaticParams.P,...
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
Sol.time_distr = toc;

if method==0
    save(strcat(folder,'/Sol',num2str(order_ct),'_PPA.mat'),'-struct','Sol');
else
    save(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'),'-struct','Sol');
end
end
end
end