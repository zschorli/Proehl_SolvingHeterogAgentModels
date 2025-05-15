% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function solves the model without aggregate risk.
%__________________________________________________________________________
function solveModel_NoAggShock(StaticParams,folder)

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
preSol = struct; 
NK = length(StaticParams.KGrid); 
NL = length(StaticParams.LGrid); 
NLt = length(StaticParams.LtaxGrid); 

%--------------------------------------------------------------------------
% Computing the solution when the aggregate state is fixed to 1 (good state)
%--------------------------------------------------------------------------     
c_g = repmat(log(StaticParams.kGrid_pol+1)+1e-4,[2,1,NK,NL,NLt]);
l_g = max(StaticParams.l_min,repmat(permute(StaticParams.LGrid,[3,4,1,2]),[2,numel(StaticParams.kGrid_pol),numel(StaticParams.KGrid),1,NLt]));
l_g(1,:,:,:,:) = StaticParams.lu;
k_prime_g = (StaticParams.wealth(1,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,...
            StaticParams.kGrid_pol,l_g)-c_g);
% c_g = StaticParams.wealth(1,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,...
%       StaticParams.kGrid_pol,l_g)-k_prime_g;

tic;
[c_g,l_g,k_prime_g,diff_g]...
= OptPol_NoAggShock_PFI(true,c_g,l_g,k_prime_g,StaticParams);
preSol.time_g = toc;

preSol.k_prime_g = k_prime_g.*(k_prime_g>eps);
preSol.c_g = c_g;
preSol.l_g = l_g;
preSol.diffs_g = diff_g;

save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');

%--------------------------------------------------------------------------
% Computing the solution when the aggregate state is fixed to 0 (bad state)
%--------------------------------------------------------------------------
c_b = repmat(log(StaticParams.kGrid_pol+1)+1e-4,[2,1,NK,NL,NLt]);
l_b = max(StaticParams.l_min,repmat(permute(StaticParams.LGrid,[3,4,1,2]),[2,numel(StaticParams.kGrid_pol),numel(StaticParams.KGrid),1,NLt]));
l_b(1,:,:,:,:) = StaticParams.lu;
k_prime_b = (StaticParams.wealth(1,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,...
            StaticParams.kGrid_pol,l_b)-c_b);
% c_b = StaticParams.wealth(1,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,...
%       StaticParams.kGrid_pol,l_b)-k_prime_b;

tic;
[c_b,l_b,k_prime_b,diff_b]...
= OptPol_NoAggShock_PFI(false,c_b,l_b,k_prime_b,StaticParams);
preSol.time_b = toc;

preSol.k_prime_b = k_prime_b.*(k_prime_b>eps);
preSol.c_b = c_b;
preSol.l_b = l_b;
preSol.diffs_b = diff_b;

save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');

clear N c k_prime y k_prime_g k_prime_b c_g c_b y_g y_b diff_g diff_b ...
      flag_g flag_b;
% preSol = load(strcat(folder,'/preSol_PFI.mat'));
  
%--------------------------------------------------------------------------
% Computing the stationary distributions for the cases without aggregate
% shocks
%--------------------------------------------------------------------------
pdf = zeros(2,2,length(StaticParams.kGrid));  
agK_ss = zeros(2,1);  
agL_ss = zeros(2,1);  
agLtaxRatio_ss = zeros(2,1);
agC_ss = zeros(2,1);    
k = zeros([2,size(preSol.k_prime_g)]);
l_aux = zeros([2,size(preSol.k_prime_g)]);
c_aux = zeros([2,size(preSol.k_prime_g)]);
k(1,:,:,:,:,:) = preSol.k_prime_b;
k(2,:,:,:,:,:) = preSol.k_prime_g;
l_aux(1,:,:,:,:,:) = preSol.l_b;
l_aux(2,:,:,:,:,:) = preSol.l_g;
c_aux(1,:,:,:,:,:) = preSol.c_b;
c_aux(2,:,:,:,:,:) = preSol.c_g;
P = zeros([2,size(StaticParams.P_g)]);
P(1,:,:) = StaticParams.P_b;
P(2,:,:) = StaticParams.P_g;   

tic;
for bool_gb=0:1 
    k_prime = max(StaticParams.k_min,squeeze(k(bool_gb+1,:,:,:,:,:)));
    l = max(StaticParams.k_min,squeeze(l_aux(bool_gb+1,:,:,:,:,:)));
    c = max(StaticParams.k_min,squeeze(c_aux(bool_gb+1,:,:,:,:,:)));

    [Yz,Yk] = ndgrid(1:2,StaticParams.kGrid_pol,1,1,1);
    if numel(StaticParams.LtaxGrid)>1
        k_prime_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,k_prime,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)),Ltax.*ones(size(Yz)));    
        l_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,l,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)),Ltax.*ones(size(Yz)));
    else
        k_prime_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,k_prime,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)));    
        l_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,l,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)));
    end
    p = StaticParams.p(2*bool_gb+(1:2))*2;
    agL_stationary = zeros(length(StaticParams.KGrid),length(StaticParams.LtaxGrid));
    for j=1:length(StaticParams.LtaxGrid)
        parfor i=1:length(StaticParams.KGrid)
        K = StaticParams.KGrid(i);
        Ltax = StaticParams.LtaxGrid(j);
        agL_stationary(i,j) = fzero(@(L) L-calcStationaryAggrL(k_prime_approx(K,L,Ltax),l_approx(K,L,Ltax),squeeze(P(bool_gb+1,:,:)),p,StaticParams),[StaticParams.LGrid(1),StaticParams.LGrid(end)]);
        end
    end
    if numel(StaticParams.LtaxGrid)>1
        k_prime_approx = @(K,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,k_prime,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,StaticParams.LtaxGrid,agL_stationary,K,Ltax).*ones(size(Yz)),Ltax.*ones(size(Yz)));
        l_approx = @(K,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,l,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,StaticParams.LtaxGrid,agL_stationary,K,Ltax).*ones(size(Yz)),Ltax.*ones(size(Yz)));
    else
        k_prime_approx = @(K,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,k_prime,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,agL_stationary,K).*ones(size(Yz)));
        l_approx = @(K,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,l,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,agL_stationary,K).*ones(size(Yz)));
    end
    agLtaxRatio_stationary = ones(size(StaticParams.KGrid));
    if numel(StaticParams.LtaxGrid)>1
        parfor i=1:length(StaticParams.KGrid)
            K = StaticParams.KGrid(i);
            agLtaxRatio_stationary(i) = fzero(@(Ltax) Ltax-calcStationaryAggrL(k_prime_approx(K,Ltax),l_approx(K,Ltax),squeeze(P(bool_gb+1,:,:)),p,StaticParams)...
                                                      ./calcStationaryAggrL(k_prime_approx(K,Ltax),l_approx(K,Ltax).^(1-StaticParams.rho),squeeze(P(bool_gb+1,:,:)),p,StaticParams),...
                                                      [StaticParams.LtaxGrid(1),StaticParams.LtaxGrid(end)]);
        end
        k_prime_approx = @(K) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,k_prime,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,StaticParams.LtaxGrid,agL_stationary,K,interp1(StaticParams.KGrid,agLtaxRatio_stationary,K)).*ones(size(Yz)),interp1(StaticParams.KGrid,agLtaxRatio_stationary,K).*ones(size(Yz)));
        l_approx = @(K) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,l,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,StaticParams.LtaxGrid,agL_stationary,K,interp1(StaticParams.KGrid,agLtaxRatio_stationary,K)).*ones(size(Yz)),interp1(StaticParams.KGrid,agLtaxRatio_stationary,K).*ones(size(Yz)));
        c_approx = @(K) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,c,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,StaticParams.LtaxGrid,agL_stationary,K,interp1(StaticParams.KGrid,agLtaxRatio_stationary,K)).*ones(size(Yz)),interp1(StaticParams.KGrid,agLtaxRatio_stationary,K).*ones(size(Yz)));
    else
        k_prime_approx = @(K) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,k_prime,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,agL_stationary,K).*ones(size(Yz)));
        l_approx = @(K) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,l,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,agL_stationary,K).*ones(size(Yz)));
        c_approx = @(K) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,c,Yz,Yk,K.*ones(size(Yz)),...
                              interpn(StaticParams.KGrid,agL_stationary,K).*ones(size(Yz)));
    end
    agK_stationary = fzero(@(K) K-calcStationaryDistr(k_prime_approx(K),l_approx(K),squeeze(P(bool_gb+1,:,:)),p,StaticParams),[StaticParams.KGrid(1),StaticParams.KGrid(end)]);
    [~,agL_stationary,pdf_stationary] = calcStationaryDistr(k_prime_approx(agK_stationary),l_approx(agK_stationary),squeeze(P(bool_gb+1,:,:)),p,StaticParams);
    [~,agLtax_stationary] = calcStationaryDistr(k_prime_approx(agK_stationary),l_approx(agK_stationary).^(1-StaticParams.rho),squeeze(P(bool_gb+1,:,:)),p,StaticParams);
    [~,agC_stationary] = calcStationaryDistr(k_prime_approx(agK_stationary),c_approx(agK_stationary),squeeze(P(bool_gb+1,:,:)),p,StaticParams);
    
    pdf(bool_gb+1,:,:) = pdf_stationary;
    agK_ss(bool_gb+1) = agK_stationary;
    agL_ss(bool_gb+1) = agL_stationary;
    agLtaxRatio_ss(bool_gb+1) = agL_stationary./agLtax_stationary;
    agC_ss(bool_gb+1) = agC_stationary;
end
preSol.time_distr = toc;

preSol.pdf_cond_b = squeeze(pdf(1,:,:)); 
preSol.pdf_cond_g = squeeze(pdf(2,:,:)); 
preSol.cdf_cond_b = cummax(min(1,cumsum(preSol.pdf_cond_b,2)),2); 
preSol.cdf_cond_g = cummax(min(1,cumsum(preSol.pdf_cond_g,2)),2);
preSol.cdf_cond_b(:,end) = 1;
preSol.cdf_cond_g(:,end) = 1;
preSol.agK_stationary_b = agK_ss(1);
preSol.agK_stationary_g = agK_ss(2);
preSol.agL_stationary_b = agL_ss(1);
preSol.agL_stationary_g = agL_ss(2);
preSol.agLtaxRatio_stationary_b = agLtaxRatio_ss(1);
preSol.agLtaxRatio_stationary_g = agLtaxRatio_ss(2);
preSol.agC_stationary_b = agC_ss(1);
preSol.agC_stationary_g = agC_ss(2);
[Yz,Yk] = ndgrid(1:2,StaticParams.kGrid,1,1,1);
if numel(StaticParams.LtaxGrid)>1
    l_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,preSol.l_b,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)),Ltax.*ones(size(Yz)));
else
    l_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,preSol.l_b,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)));
end
pdf = calcEoPDistr(1-l_approx(preSol.agK_stationary_b,preSol.agL_stationary_b,preSol.agLtaxRatio_stationary_b),preSol.cdf_cond_b,StaticParams.lGrid);
preSol.pdf_cond_l_b = pdf(:,size(pdf,2):-1:1);
if numel(StaticParams.LtaxGrid)>1
    l_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,preSol.l_g,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)),Ltax.*ones(size(Yz)));
else
    l_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,preSol.l_g,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)));
end
pdf = calcEoPDistr(1-l_approx(preSol.agK_stationary_g,preSol.agL_stationary_g,preSol.agLtaxRatio_stationary_g),preSol.cdf_cond_g,StaticParams.lGrid);
preSol.pdf_cond_l_g = pdf(:,size(pdf,2):-1:1);
if numel(StaticParams.LtaxGrid)>1
    c_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,preSol.c_b,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)),Ltax.*ones(size(Yz)));
else
    c_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,preSol.c_b,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)));
end
pdf = calcEoPDistr(c_approx(preSol.agK_stationary_b,preSol.agL_stationary_b,preSol.agLtaxRatio_stationary_b),preSol.cdf_cond_b,StaticParams.cGrid);
preSol.pdf_cond_c_b = pdf;
if numel(StaticParams.LtaxGrid)>1
    c_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,preSol.c_g,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)),Ltax.*ones(size(Yz)));
else
    c_approx = @(K,L,Ltax) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,preSol.c_g,Yz,Yk,K.*ones(size(Yz)),L.*ones(size(Yz)));
end
pdf = calcEoPDistr(c_approx(preSol.agK_stationary_g,preSol.agL_stationary_g,preSol.agLtaxRatio_stationary_g),preSol.cdf_cond_g,StaticParams.cGrid);
preSol.pdf_cond_c_g = pdf;
clear k k_prime k_prime_approx pdf pdf_cond pdf_cond_new cdf_cond agK ...
      diffs iter bool_gb i iter_ct;
  
save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');
end