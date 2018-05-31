% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This file computes the expected ergodic distributions and Euler equation 
% errors of te different algorithm solutions.
%__________________________________________________________________________
function results_statDistr(case_nr,StaticParams,folder)

% truncation order 
switch case_nr
    case 1
        M = 3;
        M_KS = 4;
    case 2
        M = 4;
        M_KS = 3;
    case 4
        M = 4;
        M_KS = 1;
    otherwise
        M = 3;
        M_KS = 1;
end

if case_nr<=2
%--------------------------------------------------------------------------
% Computing the stationary distribution for the Krusell-Smith Algo from
% Maliar, Maliar, Valli (2010)
%--------------------------------------------------------------------------
for order=1:M_KS
% Initialization___________________________________________________________
KS = load(strcat(folder,'/KS',num2str(order),'_Sol',num2str(case_nr),'.mat'));
KS.pdf_cond = zeros(4,length(StaticParams.kGrid));  
KS.pdf_cond(:,getCapIdx(StaticParams.kss,StaticParams.kGrid)) = 1;
KS.cdf_cond = cumsum(KS.pdf_cond,2);
KS.cdf_cond(:,end) = 1;
KS.pdf = [[StaticParams.ur(1),StaticParams.er(1)]*KS.pdf_cond(1:2,:);...
          [StaticParams.ur(2),StaticParams.er(2)]*KS.pdf_cond(3:4,:)];
diffs_p = 1e8*ones(2,1);
iter = 1;   

% Iteration________________________________________________________________
while max(diffs_p(end)) > StaticParams.criter_pdf
    K = KS.pdf*StaticParams.kGrid';
    switch order
        case 1
        k_prime = max(0,[interpn(KS.k,KS.km,KS.kprime(:,:,1,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))))';...
                     interpn(KS.k,KS.km,KS.kprime(:,:,1,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))))';...
                     interpn(KS.k,KS.km,KS.kprime(:,:,2,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))))';...
                     interpn(KS.k,KS.km,KS.kprime(:,:,2,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))))']);
        case 2
        sd = [sqrt(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^2);sqrt(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^2)];
        k_prime = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,1,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))))';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,1,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))))';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,2,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))))';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,2,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))))']);
        case 3
        sd = [sqrt(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^2);sqrt(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^2)];
        skew = [(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^3)/sd(1)^3;(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^3)/sd(2)^3];
        k_prime = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,1,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),min(KS.kskew_max,max(KS.kskew_min,skew(1))))';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,1,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),min(KS.kskew_max,max(KS.kskew_min,skew(1))))';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,2,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),min(KS.kskew_max,max(KS.kskew_min,skew(2))))';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,2,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),min(KS.kskew_max,max(KS.kskew_min,skew(2))))']);
        case 4
        sd = [sqrt(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^2);sqrt(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^2)];
        skew = [(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^3)/sd(1)^3;(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^3)/sd(2)^3];
        kurt = [(KS.pdf(1,:)*(StaticParams.kGrid'-K(1)).^4)/sd(1)^4;(KS.pdf(2,:)*(StaticParams.kGrid'-K(2)).^4)/sd(2)^4];
        k_prime = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,1,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),min(KS.kskew_max,max(KS.kskew_min,skew(1))),min(KS.kkurt_max,max(KS.kkurt_min,kurt(1))))';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,1,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(1))),min(KS.kvar_max,max(KS.kvar_min,sd(1))),min(KS.kskew_max,max(KS.kskew_min,skew(1))),min(KS.kkurt_max,max(KS.kkurt_min,kurt(1))))';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,2,1),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),min(KS.kskew_max,max(KS.kskew_min,skew(2))),min(KS.kkurt_max,max(KS.kkurt_min,kurt(2))))';...
                     interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,2,2),StaticParams.kGrid_pol,min(KS.km_max,max(KS.km_min,K(2))),min(KS.kvar_max,max(KS.kvar_min,sd(2))),min(KS.kskew_max,max(KS.kskew_min,skew(2))),min(KS.kkurt_max,max(KS.kkurt_min,kurt(2))))']);
    end                 
    pdf_cond_new = calcEoPDistr(interp1(StaticParams.kGrid_pol,k_prime',StaticParams.kGrid)',...
                   KS.cdf_cond,StaticParams.kGrid);
    pdf_cond_new = (repmat(1./(StaticParams.p'*StaticParams.P)',[1,4])...
                   .*StaticParams.P')...
                   *(repmat(StaticParams.p,[1,size(KS.pdf_cond,2)])...
                   .*pdf_cond_new);
    pdf_cond_new = pdf_cond_new.*(pdf_cond_new>eps);
    pdf_cond = pdf_cond_new./repmat(sum(pdf_cond_new,2),[1,size(pdf_cond_new,2)]);
    pdf = [[StaticParams.ur(1),StaticParams.er(1)]*pdf_cond_new(1:2,:);...
           [StaticParams.ur(2),StaticParams.er(2)]*pdf_cond_new(3:4,:)];
    diffs_p(iter) = norm(pdf_cond_new-KS.pdf_cond,2);
    KS.pdf_cond = pdf_cond_new;
    KS.pdf = pdf;
    KS.cdf_cond = cummax(min(1,cumsum(pdf_cond,2)),2);
    KS.cdf_cond(1,KS.cdf_cond(1,:)==max(KS.cdf_cond(1,:))) = 1;
    KS.cdf_cond(2,KS.cdf_cond(2,:)==max(KS.cdf_cond(2,:))) = 1;
    iter = iter+1;
end
KS.cdf = cumsum(KS.pdf,2);
KS.cdf(:,end) = 1;
clear k_prime c_prime pdf_cond_new cdf_new weights diffs_p iter ;
save(strcat(folder,'/KS',num2str(order),'_Sol',num2str(case_nr),'.mat'),'-struct','KS');
end


%--------------------------------------------------------------------------
% Converting the endogenous grid from the Reiter (backward induction) Algo
% from Reiter (2010) to our fixed grid for comparison
%--------------------------------------------------------------------------
R = load(strcat(folder,'/R_Sol',num2str(case_nr),'.mat'));
wealthR = repmat(R.VAac.x',[R.Params2.nGridMom,1,2,2])+1./R.VAac.dv;
rateR = StaticParams.alpha.*...
       repmat(permute([1-StaticParams.delta_a;1+StaticParams.delta_a],[3,4,1,2]),...
              [R.Params2.nGridMom,R.Params2.n,1,2])...
       .*(repmat(R.Params2.GridMom(:,1),[1,R.Params2.n,2,2])...
          ./(StaticParams.l_bar.*...
             repmat(permute(StaticParams.er,[3,4,1,2]),...
                    [R.Params2.nGridMom,R.Params2.n,1,2])))...
       .^(StaticParams.alpha-1);  
wageR = (1-StaticParams.alpha).*...
       repmat(permute([1-StaticParams.delta_a;1+StaticParams.delta_a],[3,4,1,2]),...
              [R.Params2.nGridMom,R.Params2.n,1,2])...
       .*(repmat(R.Params2.GridMom(:,1),[1,R.Params2.n,2,2])...
          ./(StaticParams.l_bar.*...
             repmat(permute(StaticParams.er,[3,4,1,2]),...
                    [R.Params2.nGridMom,R.Params2.n,1,2])))...
       .^(StaticParams.alpha);  
% compute the endogenous start capital   
kGrid = (wealthR - ...
        ((StaticParams.l_bar-StaticParams.mu.*...
          repmat(permute(StaticParams.ur./StaticParams.er,[3,4,1,2]),...
                 [R.Params2.nGridMom,R.Params2.n,1,2]))...
         .*repmat(permute([0,1],[3,4,1,2]),[R.Params2.nGridMom,R.Params2.n,2,1])...
        +StaticParams.mu.*repmat(permute([1,0],[3,4,1,2]),[R.Params2.nGridMom,R.Params2.n,2,1]))...
        .*wageR)...
        ./(1-StaticParams.delta+rateR);
kGrid(:,2:end,:,1) = kGrid(:,1:end-1,:,1);    
kGrid(:,1,:,1) = 0;
% compute the corresponding policies
c_prime_endGrid = 1./R.VAac.dv;
c_prime_endGrid(:,2:end,:,1) = c_prime_endGrid(:,1:end-1,:,1);
c_prime_endGrid(:,1,:,1) = StaticParams.mu.*wageR(:,1,:,1);
k_prime_endGrid = repmat(R.VAac.x',[R.Params2.nGridMom,1,2,2]);
k_prime_endGrid(:,2:end,:,1) = k_prime_endGrid(:,1:end-1,:,1);    
k_prime_endGrid(:,1,:,1) = 0;
% convert grid to comparison grid
R.c_prime = zeros(R.Params2.nGridMom,4,length(StaticParams.kGrid_pol));
R.k_prime = zeros(R.Params2.nGridMom,4,length(StaticParams.kGrid_pol));
R.idx = (1:R.Params2.nGridMom)';
for i=1:R.Params2.nGridMom
for i_ag = 1:2
for i_id = 1:2
R.c_prime(i,2*(i_ag-1)+i_id,:) = interp1(squeeze(kGrid(i,:,i_ag,i_id)),...
    squeeze(c_prime_endGrid(i,:,i_ag,i_id)),StaticParams.kGrid_pol,'linear','extrap');
R.k_prime(i,2*(i_ag-1)+i_id,:) = interp1(squeeze(kGrid(i,:,i_ag,i_id)),...
    squeeze(k_prime_endGrid(i,:,i_ag,i_id)),StaticParams.kGrid_pol,'linear','extrap');
end
end
end
clear wealthR rateR wageR kGrid c_prime_endGrid k_prime_endGrid i j i_ag i_id;

%--------------------------------------------------------------------------
% Computing the stationary distribution for the Reiter Algo (backward 
% induction) from Reiter (2010)
%--------------------------------------------------------------------------
% Initialization___________________________________________________________
R.pdf_cond = zeros(4,length(StaticParams.kGrid));  
R.pdf_cond(:,getCapIdx(40,StaticParams.kGrid)) = 1;
R.cdf_cond = cumsum(R.pdf_cond,2);
R.cdf_cond(:,end) = 1;
R.pdf = [[StaticParams.ur(1),StaticParams.er(1)]*R.pdf_cond(1:2,:);...
          [StaticParams.ur(2),StaticParams.er(2)]*R.pdf_cond(3:4,:)];
K = R.pdf*StaticParams.kGrid';
diffs_p = 1e8*ones(2,1);
iter = 1;   

% Iteration________________________________________________________________
while max(diffs_p(end)) > StaticParams.criter_pdf
    k_prime = max(0,[getCApprox(R.Params2.GridMat,squeeze(R.k_prime(:,1,:)),R.idx,K(1));...
                     getCApprox(R.Params2.GridMat,squeeze(R.k_prime(:,2,:)),R.idx,K(1));...
                     getCApprox(R.Params2.GridMat,squeeze(R.k_prime(:,3,:)),R.idx,K(2));...
                     getCApprox(R.Params2.GridMat,squeeze(R.k_prime(:,4,:)),R.idx,K(2))]);
    pdf_cond_new = calcEoPDistr(interp1(StaticParams.kGrid_pol,k_prime',StaticParams.kGrid)',...
                   R.cdf_cond,StaticParams.kGrid);
    pdf_cond_new = (repmat(1./(StaticParams.p'*StaticParams.P)',[1,4])...
                   .*StaticParams.P')...
                   *(repmat(StaticParams.p,[1,size(R.pdf_cond,2)])...
                   .*pdf_cond_new);
    pdf_cond_new = pdf_cond_new.*(pdf_cond_new>eps);
    pdf_cond = pdf_cond_new./repmat(sum(pdf_cond_new,2),[1,size(pdf_cond_new,2)]);
    pdf = [[StaticParams.ur(1),StaticParams.er(1)]*pdf_cond_new(1:2,:);...
           [StaticParams.ur(2),StaticParams.er(2)]*pdf_cond_new(3:4,:)];
    K = pdf*StaticParams.kGrid';       
    diffs_p(iter) = norm(pdf_cond_new-R.pdf_cond,2);
    R.pdf_cond = pdf_cond_new;
    R.pdf = pdf;
    R.cdf_cond = cummax(min(1,cumsum(pdf_cond,2)),2);
    R.cdf_cond(1,R.cdf_cond(1,:)==max(R.cdf_cond(1,:))) = 1;
    R.cdf_cond(2,R.cdf_cond(2,:)==max(R.cdf_cond(2,:))) = 1;
    iter = iter+1;
end
R.cdf = cumsum(R.pdf,2);
R.cdf(:,end) = 1;
clear k_prime c_prime pdf_cond_new cdf_new weights diffs_p iter ;
save(strcat(folder,'/R_Sol',num2str(case_nr),'.mat'),'-struct','R');


%--------------------------------------------------------------------------
% Computing the stationary distribution for the Den Haan-Rendahl Algo from
% Den Haan and Rendahl (2010)
%--------------------------------------------------------------------------
% Initialization___________________________________________________________
DR = load(strcat(folder,'/DR_Sol',num2str(case_nr),'.mat'));
DR.pdf_cond = zeros(4,length(StaticParams.kGrid));  
DR.pdf_cond(:,getCapIdx(StaticParams.kss,StaticParams.kGrid)) = 1;
DR.cdf_cond = cumsum(DR.pdf_cond,2);
DR.cdf_cond(:,end) = 1;
DR.pdf = [[StaticParams.ur(1),StaticParams.er(1)]*DR.pdf_cond(1:2,:);...
          [StaticParams.ur(2),StaticParams.er(2)]*DR.pdf_cond(3:4,:)];
diffs_p = 1e8*ones(2,1);
iter = 1;   

% Iteration________________________________________________________________
while max(diffs_p(end)) > StaticParams.criter_pdf
    K = DR.pdf_cond*StaticParams.kGrid';
    k_prime = max(0,[interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,2,2),StaticParams.kGrid_pol,max(min(K(2),DR.Kemax),DR.Kemin),max(min(K(1),DR.Kumax),DR.Kumin))';...
                     interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,1,2),StaticParams.kGrid_pol,max(min(K(2),DR.Kemax),DR.Kemin),max(min(K(1),DR.Kumax),DR.Kumin))';...
                     interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,2,1),StaticParams.kGrid_pol,max(min(K(4),DR.Kemax),DR.Kemin),max(min(K(3),DR.Kumax),DR.Kumin))';...
                     interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,1,1),StaticParams.kGrid_pol,max(min(K(4),DR.Kemax),DR.Kemin),max(min(K(3),DR.Kumax),DR.Kumin))']);
    pdf_cond_new = calcEoPDistr(interp1(StaticParams.kGrid_pol,k_prime',StaticParams.kGrid)',...
                   DR.cdf_cond,StaticParams.kGrid);
    pdf_cond_new = (repmat(1./(StaticParams.p'*StaticParams.P)',[1,4])...
                   .*StaticParams.P')...
                   *(repmat(StaticParams.p,[1,size(DR.pdf_cond,2)])...
                   .*pdf_cond_new);
    pdf_cond_new = pdf_cond_new.*(pdf_cond_new>eps);
    pdf_cond = pdf_cond_new./repmat(sum(pdf_cond_new,2),[1,size(pdf_cond_new,2)]);
    pdf = [[StaticParams.ur(1),StaticParams.er(1)]*pdf_cond_new(1:2,:);...
           [StaticParams.ur(2),StaticParams.er(2)]*pdf_cond_new(3:4,:)];
    diffs_p(iter) = norm(pdf_cond_new-DR.pdf_cond,2);
    DR.pdf_cond = pdf_cond_new;
    DR.pdf = pdf;
    DR.cdf_cond = cummax(min(1,cumsum(pdf_cond,2)),2);
    DR.cdf_cond(1,DR.cdf_cond(1,:)==max(DR.cdf_cond(1,:))) = 1;
    DR.cdf_cond(2,DR.cdf_cond(2,:)==max(DR.cdf_cond(2,:))) = 1;
    iter = iter+1;
end
DR.cdf = cumsum(DR.pdf,2);
DR.cdf(:,end) = 1;
clear k_prime c_prime pdf_cond_new cdf_new weights diffs_p iter ;
save(strcat(folder,'/DR_Sol',num2str(case_nr),'.mat'),'-struct','DR');
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

diffs_p = 1e8*ones(2,1);
iter = 1;
w_bounds = repmat([1e8;-1e8],[1,Poly.phi_N+1]);


% Iteration________________________________________________________________
tic;
while max(diffs_p(end)) > StaticParams.criter_pdf && iter<10001
    k_prime = max(0,[getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,1,:)),Sol.idx,weights(1,:));...
                     getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2,:)),Sol.idx,weights(1,:));...
                     getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,3,:)),Sol.idx,weights(2,:));...
                     getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,4,:)),Sol.idx,weights(2,:))]);
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
    if iter==10001
       Sol.DistrNoConvergence = {true,diffs_p(iter-1)}; 
    end
end
Sol.pdf = [[StaticParams.ur(1),StaticParams.er(1)]*Sol.pdf_cond(1:2,:);...
           [StaticParams.ur(2),StaticParams.er(2)]*Sol.pdf_cond(3:4,:)];
Sol.cdf_cond = cdf_new;
Sol.cdf = cummax(min(1,cumsum(Sol.pdf,2)),2);
Sol.cdf(1,Sol.cdf(1,:)==max(Sol.cdf(1,:))) = 1;
Sol.cdf(2,Sol.cdf(2,:)==max(Sol.cdf(2,:))) = 1;
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