% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This function solves the model without aggregate risk.
%__________________________________________________________________________
function solveModel_NoAggShock(StaticParams,folder)

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
preSol = struct;

%--------------------------------------------------------------------------
% Computing the solution when the aggregate state is fixed to 1 (good state)
%--------------------------------------------------------------------------     
F_h = @(c) min((StaticParams.eta.*StaticParams.uq_large.^(-StaticParams.eps).*c)...
          ,(1/(1-StaticParams.theta).*repmat(repmat(StaticParams.aGrid_pol,[2,1]),[1,1,StaticParams.Np,StaticParams.NA])...
          ./StaticParams.q_large));

c_g = 0.1.*StaticParams.wealth(1);
h_g = F_h(c_g);
a_prime_g = StaticParams.wealth(1)-c_g-StaticParams.uq_large.*h_g;

tstart = cputime;
[c_g,h_g,a_prime_g,diff_g]...
= OptPol_NoAggShock_PFI(true,c_g,h_g,a_prime_g,StaticParams);
preSol.time_g = cputime-tstart;

preSol.a_prime_g = a_prime_g;
preSol.h_g = h_g;
preSol.c_g = c_g;
preSol.diffs_g = diff_g;

save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');

%--------------------------------------------------------------------------
% Computing the solution when the aggregate state is fixed to 0 (bad state)
%--------------------------------------------------------------------------
c_b = 0.1.*StaticParams.wealth(0);
h_b = F_h(c_b);
a_prime_b = StaticParams.wealth(0)-c_b-StaticParams.uq_large.*h_b;

tstart = cputime;
[c_b,h_b,a_prime_b,diff_b]...
= OptPol_NoAggShock_PFI(false,c_b,h_b,a_prime_b,StaticParams);
preSol.time_b = cputime-tstart;

preSol.a_prime_b = a_prime_b;
preSol.c_b = c_b;
preSol.h_b = h_b;
preSol.diffs_b = diff_b;

save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');

clear N c a_prime y a_prime_g a_prime_b c_g c_b y_g y_b diff_g diff_b ...
      flag_g flag_b;
  
%--------------------------------------------------------------------------
% Computing the stationary distributions for the cases without aggregate
% shocks
%--------------------------------------------------------------------------
pdf = zeros(2,2,length(StaticParams.aGrid));  
agA_ss = zeros(2,1);  
agR_ss = zeros(2,1);  
agC_ss = zeros(2,1);    
a = zeros([2,size(preSol.a_prime_g)]);
h_aux = zeros([2,size(preSol.a_prime_g)]);
c_aux = zeros([2,size(preSol.a_prime_g)]);
a(1,:,:,:,:,:) = preSol.a_prime_b;
a(2,:,:,:,:,:) = preSol.a_prime_g;
h_aux(1,:,:,:,:,:) = preSol.h_b;
h_aux(2,:,:,:,:,:) = preSol.h_g;
c_aux(1,:,:,:,:,:) = preSol.c_b;
c_aux(2,:,:,:,:,:) = preSol.c_g;
P = zeros([2,size(StaticParams.P_g)]);
P(1,:,:) = StaticParams.P_b;
P(2,:,:) = StaticParams.P_g; 


tstart = cputime;
for bool_gb=0:1 
    a_prime = max(StaticParams.a_min,squeeze(a(bool_gb+1,:,:,:,:,:)));
    h = max(0,squeeze(h_aux(bool_gb+1,:,:,:,:,:)));
    c = max(0,squeeze(c_aux(bool_gb+1,:,:,:,:,:)));

    [Yz,Ya] = ndgrid(1:2,StaticParams.aGrid_pol,1,1,1);
    a_prime_approx = @(R,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,a_prime,Yz,Ya,R.*ones(size(Yz)),q.*ones(size(Yz)));    
    h_approx = @(R,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,h,Yz,Ya,R.*ones(size(Yz)),q.*ones(size(Yz)));
    c_approx = @(R,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,c,Yz,Ya,R.*ones(size(Yz)),q.*ones(size(Yz)));
    
    p = StaticParams.p(2*bool_gb+(1:2))*2;

    agR_stationary = zeros(length(StaticParams.qGrid),1);
    marketCl = zeros(length(StaticParams.qGrid),1);
    parfor i=1:length(StaticParams.qGrid)
        pr = StaticParams.qGrid(i);
        agR_stationary(i) = fzero(@(R) pr-calcStationaryDistr(a_prime_approx(R,pr),h_approx(R,pr),squeeze(P(bool_gb+1,:,:)),p,StaticParams),[StaticParams.RGrid(1),StaticParams.RGrid(end)]);
        marketCl(i) = 1-calcStationaryAggrH(a_prime_approx(agR_stationary(i),pr),h_approx(agR_stationary(i),pr),squeeze(P(bool_gb+1,:,:)),p,StaticParams);
    end
    agPr_stationary = fzero(@(q) interpn(StaticParams.qGrid,marketCl,q,'pchip'),[StaticParams.qGrid(1),StaticParams.qGrid(end)]);
    agR_stationary = fzero(@(R) agPr_stationary-calcStationaryDistr(a_prime_approx(R,agPr_stationary),h_approx(R,agPr_stationary),squeeze(P(bool_gb+1,:,:)),p,StaticParams),[StaticParams.RGrid(1),StaticParams.RGrid(end)]);
    
    [~,~,pdf_stationary] = calcStationaryDistr(a_prime_approx(agR_stationary,agPr_stationary),h_approx(agR_stationary,agPr_stationary),squeeze(P(bool_gb+1,:,:)),p,StaticParams);
    [~,agC_stationary] = calcStationaryDistr(a_prime_approx(agR_stationary,agPr_stationary),c_approx(agR_stationary,agPr_stationary),squeeze(P(bool_gb+1,:,:)),p,StaticParams);
    
    pdf(bool_gb+1,:,:) = pdf_stationary;
    agA_ss(bool_gb+1) = agPr_stationary;
    agR_ss(bool_gb+1) = agR_stationary-1;
    agC_ss(bool_gb+1) = agC_stationary;
end

preSol.pdf_cond_b = squeeze(pdf(1,:,:)); 
preSol.pdf_cond_g = squeeze(pdf(2,:,:)); 
preSol.cdf_cond_b = cummax(min(1,cumsum(preSol.pdf_cond_b,2)),2); 
preSol.cdf_cond_g = cummax(min(1,cumsum(preSol.pdf_cond_g,2)),2);
preSol.cdf_cond_b(:,end) = 1;
preSol.cdf_cond_g(:,end) = 1;

preSol.pdf_b = [StaticParams.ur(1),StaticParams.er(1)]*preSol.pdf_cond_b;
preSol.cdf_b = cummax(min(1,cumsum(preSol.pdf_b,2)),2);
preSol.cdf_b(1,preSol.cdf_b==max(preSol.cdf_b)) = 1;
preSol.Lorenz_b=(tril(repmat(preSol.pdf_b,[size(preSol.pdf_b,2),1]))*StaticParams.aGrid')'/(preSol.pdf_b*StaticParams.aGrid');
preSol.gini_b = sum(sum((preSol.pdf_b'*preSol.pdf_b).*abs(repmat(StaticParams.aGrid',[1,size(preSol.pdf_b,2)])-repmat(StaticParams.aGrid,[size(preSol.pdf_b,2),1])),1),2)...
                /(2*StaticParams.aGrid*preSol.pdf_b');
preSol.pdf_g = [StaticParams.ur(1),StaticParams.er(1)]*preSol.pdf_cond_g;
preSol.cdf_g = cummax(min(1,cumsum(preSol.pdf_g,2)),2);
preSol.cdf_g(1,preSol.cdf_g==max(preSol.cdf_g)) = 1;
preSol.Lorenz_g=(tril(repmat(preSol.pdf_g,[size(preSol.pdf_g,2),1]))*StaticParams.aGrid')'/(preSol.pdf_g*StaticParams.aGrid');
preSol.gini_g = sum(sum((preSol.pdf_g'*preSol.pdf_g).*abs(repmat(StaticParams.aGrid',[1,size(preSol.pdf_g,2)])-repmat(StaticParams.aGrid,[size(preSol.pdf_g,2),1])),1),2)...
                /(2*StaticParams.aGrid*preSol.pdf_g');

preSol.agA_stationary_b = agA_ss(1);
preSol.agA_stationary_g = agA_ss(2);
preSol.agR_stationary_b = agR_ss(1);
preSol.agR_stationary_g = agR_ss(2);
preSol.agC_stationary_b = agC_ss(1);
preSol.agC_stationary_g = agC_ss(2);

[Yz,Ya] = ndgrid(1:2,StaticParams.aGrid,1,1,1);
c_approx = @(r,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,preSol.c_b,Yz,Ya,(1+r).*ones(size(Yz)),q.*ones(size(Yz)));
pdf = calcEoPDistr(c_approx(preSol.agR_stationary_b,preSol.agA_stationary_b),preSol.cdf_cond_b,StaticParams.cGrid);
preSol.pdf_cond_c_b = pdf;
h_approx = @(r,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,preSol.h_b,Yz,Ya,(1+r).*ones(size(Yz)),q.*ones(size(Yz)));
pdf = calcEoPDistr(h_approx(preSol.agR_stationary_b,preSol.agA_stationary_b),preSol.cdf_cond_b,StaticParams.hGrid);
preSol.pdf_cond_h_b = pdf;
a_prime_approx = @(r,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,preSol.a_prime_b,Yz,Ya,(1+r).*ones(size(Yz)),q.*ones(size(Yz)));
b_approx = @(r,q) a_prime_approx(r,q)-q.*h_approx(r,q);
pdf = calcEoPDistr(b_approx(preSol.agR_stationary_b,preSol.agA_stationary_b),preSol.cdf_cond_b,StaticParams.bGrid);
preSol.pdf_cond_b_b = pdf;

c_approx = @(r,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,preSol.c_g,Yz,Ya,(1+r).*ones(size(Yz)),q.*ones(size(Yz)));
pdf = calcEoPDistr(c_approx(preSol.agR_stationary_g,preSol.agA_stationary_g),preSol.cdf_cond_g,StaticParams.cGrid);
preSol.pdf_cond_c_g = pdf;
h_approx = @(r,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,preSol.h_g,Yz,Ya,(1+r).*ones(size(Yz)),q.*ones(size(Yz)));
pdf = calcEoPDistr(h_approx(preSol.agR_stationary_g,preSol.agA_stationary_g),preSol.cdf_cond_g,StaticParams.hGrid);
preSol.pdf_cond_h_g = pdf;
a_prime_approx = @(r,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,preSol.a_prime_g,Yz,Ya,(1+r).*ones(size(Yz)),q.*ones(size(Yz)));
b_approx = @(r,q) a_prime_approx(r,q)-q.*h_approx(r,q);
pdf = calcEoPDistr(b_approx(preSol.agR_stationary_g,preSol.agA_stationary_g),preSol.cdf_cond_g,StaticParams.bGrid);
preSol.pdf_cond_b_g = pdf;

preSol.time_distr = cputime-tstart;  

save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');
end