% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This file computes the expected ergodic distributions of the different 
% algorithm solutions.
%__________________________________________________________________________
function results_Errors(case_nr,StaticParams,folder)

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

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
rate = @(z_ag,K) StaticParams.alpha.*...
       (1-(1-z_ag)*StaticParams.delta_a+z_ag*StaticParams.delta_a)*ones(2,1)...
       .*(repmat(K,[2,1])...
          ./(StaticParams.l_bar.*StaticParams.er(z_ag+1)*ones(2,1)))...
       .^(StaticParams.alpha-1);   
wage = @(z_ag,K) (1-StaticParams.alpha).*...
       (1-(1-z_ag)*StaticParams.delta_a+z_ag*StaticParams.delta_a)*ones(2,1)...
       .*(repmat(K,[2,1])...
          ./(StaticParams.l_bar.*StaticParams.er(z_ag+1)*ones(2,1)))...
       .^(StaticParams.alpha);   
wealth = @(z_ag,K,grid) ((StaticParams.l_bar-...
          (StaticParams.mu.*StaticParams.ur(z_ag+1)/StaticParams.er(z_ag+1)))...
          *repmat([0;1],[1,length(grid)])...
          +StaticParams.mu.*repmat([1;0],[1,length(grid)]))...
         .*repmat(wage(z_ag,K),[1,length(grid)])...
         +(1-StaticParams.delta+repmat(rate(z_ag,K),[1,length(grid)]))...
          .*repmat(grid,[2,1]);
      
R = load(strcat(folder,'/R_Sol',num2str(case_nr),'.mat'));
DR = load(strcat(folder,'/DR_Sol',num2str(case_nr),'.mat'));


% Simulate the aggregate shock_____________________________________________
% Keep this section commented out to reproduce the results from the paper!
% 
% Res = struct;
% Res.nsamples = 3000; 
% pdf = @(z) 0.5;
% proppdf = @(z_prime,z) sum(StaticParams.P(2*(z>0)+1,2*(z_prime>0)+(1:2)),2);
% proprnd = @(z) (rand>=sum(StaticParams.P(2*(z>0)+1,1:2),2));
% Res.z_ag = mhsample(0,Res.nsamples,'pdf',pdf,'proppdf',proppdf, 'proprnd',proprnd);
% clear proppdf proprnd pdf
% save('res_case1/Res.mat','-struct','Res');

Res = load(strcat(folder,'/Res.mat'));
Res.kGrid = unique([StaticParams.kGrid_pol,linspace(StaticParams.k_min,StaticParams.k_max,1801)]);
algo_idx = logical(sum(repmat(Res.kGrid,[length(StaticParams.kGrid_pol),1])==repmat(StaticParams.kGrid_pol',[1,length(Res.kGrid)]),1));
preSol = load(strcat(folder,'/preSol_PPA.mat'));
pdf_start = [[0.5,0.5]*[preSol.pdf_cond_b(1,:);preSol.pdf_cond_g(1,:)];...
             [0.5,0.5]*[preSol.pdf_cond_b(2,:);preSol.pdf_cond_g(2,:)]];
cdf_start = cummax(min(1,cumsum(pdf_start,2)),2);
cdf_start(1,cdf_start(1,:)==max(cdf_start(1,:)))=1;
cdf_start(2,cdf_start(2,:)==max(cdf_start(2,:)))=1;

v = VideoWriter(strcat(folder,'/Distr_Sim.avi'));
open(v);
k_max=120;
x=[0 0 k_max k_max];
y=[-5e-5 0 0 -5e-5];
figure
p4=fill(x,y,'g','EdgeColor','none');
set(gca,'nextplot','add');
p5=fill(x,y,'r','EdgeColor','none');
set(gca,'nextplot','add');
p1=plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*pdf_start,'m','LineWidth',2);
set(gca,'nextplot','add');
p2=plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*pdf_start,'c','LineWidth',2);
set(gca,'nextplot','add');
p3=plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*pdf_start,'b','LineWidth',2);
set(gca,'nextplot','replacechildren');
legend([p1,p2,p3,p4,p5],'Krusell-Smith', 'Reiter', 'Proximal Point Algo', 'boom', 'recession','Location','best')
legend('boxoff')
axis([0 k_max -5e-5 3.0025e-3])
title('Simulation of Capital Distributions')
xlabel('k')
ylabel('probability')
frame = getframe(gcf);
writeVideo(v,frame);

% Set variables for Krusell-Smith       
pdf_KS = repmat(permute(pdf_start,[3,1,2]),[M_KS,1,1]);
cdf_KS = repmat(permute(cdf_start,[3,1,2]),[M_KS,1,1]);
dyn_pdf_KS = pdf_KS;
dyn_cdf_KS = cdf_KS;
K_KS = zeros(M_KS,Res.nsamples,1);
sd_KS = zeros(M_KS,Res.nsamples,1);
skew_KS = zeros(M_KS,Res.nsamples,1);
kurt_KS = zeros(M_KS,Res.nsamples,1);
gini_KS = zeros(M_KS,Res.nsamples,1);
dyn_euler_K_KS = zeros(M_KS,Res.nsamples,1);
K_pr_LoM_KS = zeros(M_KS,Res.nsamples,1);
idx_posProb_KS = zeros(M_KS,2,length(Res.kGrid));
idx = max(find(pdf_KS(1,1,:)>eps,1,'last'),find(pdf_KS(1,2,:)>eps,1,'last'));
idx_posProb_KS(:,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
dyn_idx_posProb_KS = idx_posProb_KS;

% Set variables for Reiter
pdf_R = pdf_start;
cdf_R = cdf_start;
dyn_pdf_R = pdf_R;
dyn_cdf_R = cdf_R;
K_R = zeros(Res.nsamples,1);
sd_R = zeros(Res.nsamples,1);
skew_R = zeros(Res.nsamples,1);
kurt_R = zeros(Res.nsamples,1);
gini_R = zeros(Res.nsamples,1);
dyn_euler_K_R = zeros(Res.nsamples,1);
idx_posProb_R = squeeze(idx_posProb_KS(1,:,:));
dyn_idx_posProb_R = idx_posProb_R;

% Set variables for Den Haan-Rendahl
pdf_DR = pdf_start;
cdf_DR = cdf_start;
dyn_pdf_DR = pdf_DR;
dyn_cdf_DR = cdf_DR;
K_DR = zeros(Res.nsamples,1);
sd_DR = zeros(Res.nsamples,1);
skew_DR = zeros(Res.nsamples,1);
kurt_DR = zeros(Res.nsamples,1);
gini_DR = zeros(Res.nsamples,1);
K_pr_LoM_DR = zeros(Res.nsamples,1);
dyn_euler_K_DR = zeros(Res.nsamples,1);
idx_posProb_DR = idx_posProb_R;
dyn_idx_posProb_DR = idx_posProb_R;

% Set variables for PPA order 0 to M
pdf = repmat(permute(pdf_start,[4,3,1,2]),[2,M+1,1,1]);
cdf = repmat(permute(cdf_start,[4,3,1,2]),[2,M+1,1,1]);
weights = zeros(2,M+1,M+1);
for method = 0:1
for order_ct=0:M
    if method==0
        Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PPA.mat'));
    else
        Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
    end
    weights(method+1,order_ct+1,1:(Poly.phi_N+1))...
    = calcWeights(squeeze(cdf(method+1,order_ct+1,:,:)),Poly,StaticParams,0);
end
end
dyn_pdf = pdf;
dyn_cdf = cdf;
dyn_weights = weights;
K = zeros(2,M+1,Res.nsamples,1);
sd = zeros(2,M+1,Res.nsamples,1);
skew = zeros(2,M+1,Res.nsamples,1);
kurt = zeros(2,M+1,Res.nsamples,1);
gini = zeros(2,M+1,Res.nsamples,1);
K_pr_LoM = zeros(2,M+1,Res.nsamples,1);
dyn_euler_K = zeros(2,M+1,Res.nsamples,1);
idx_posProb = repmat(permute(idx_posProb_R,[3,4,1,2]),[2,M+1]);
dyn_idx_posProb = idx_posProb;

Res.stErr_KS1 = [];
if case_nr<3
    Res.stErr_KS2 = [];
    Res.stErr_KS3 = [];
end
if case_nr == 1
    Res.stErr_KS4 = [];
end
Res.stErr_R = [];
Res.stErr_DR = [];
Res.stErr_PPA0 = [];
Res.stErr_PPA1 = [];
Res.stErr_PPA2 = [];
Res.stErr_PPA3 = [];
if case_nr==2 || case_nr>3
   Res.stErr_PPA4 = [];
end
Res.stErr_PFI0 = [];
Res.stErr_PFI1 = [];
Res.stErr_PFI2 = [];
Res.stErr_PFI3 = [];
if case_nr==2 || case_nr>3
   Res.stErr_PFI4 = [];
end

Res.dynErr_KS1 = [];
if case_nr<3
    Res.dynErr_KS2 = [];
    Res.dynErr_KS3 = [];
end
if case_nr == 1
    Res.dynErr_KS4 = [];
end
Res.dynErr_R = [];
Res.dynErr_DR = [];
Res.dynErr_PPA0 = [];
Res.dynErr_PPA1 = [];
Res.dynErr_PPA2 = [];
Res.dynErr_PPA3 = [];
if case_nr==2 || case_nr>3
   Res.dynErr_PPA4 = [];
end
Res.dynErr_PFI0 = [];
Res.dynErr_PFI1 = [];
Res.dynErr_PFI2 = [];
Res.dynErr_PFI3 = [];
if case_nr==2 || case_nr>3
   Res.dynErr_PFI4 = [];
end

Res.stErr_KS1_algoGrid = [];
if case_nr<3
    Res.stErr_KS2_algoGrid = [];
    Res.stErr_KS3_algoGrid = [];
end
if case_nr == 1
    Res.stErr_KS4_algoGrid = [];
end
Res.stErr_R_algoGrid = [];
Res.stErr_DR_algoGrid = [];
Res.stErr_PPA0_algoGrid = [];
Res.stErr_PPA1_algoGrid = [];
Res.stErr_PPA2_algoGrid = [];
Res.stErr_PPA3_algoGrid = [];
if case_nr==2 || case_nr>3
   Res.stErr_PPA4_algoGrid = [];
end
Res.stErr_PFI0_algoGrid = [];
Res.stErr_PFI1_algoGrid = [];
Res.stErr_PFI2_algoGrid = [];
Res.stErr_PFI3_algoGrid = [];
if case_nr==2 || case_nr>3
   Res.stErr_PFI4_algoGrid = [];
end

Res.dynErr_KS1_algoGrid = [];
if case_nr<3
    Res.dynErr_KS2_algoGrid = [];
    Res.dynErr_KS3_algoGrid = [];
end
if case_nr == 1
    Res.dynErr_KS4_algoGrid = [];
end
Res.dynErr_R_algoGrid = [];
Res.dynErr_DR_algoGrid = [];
Res.dynErr_PPA0_algoGrid = [];
Res.dynErr_PPA1_algoGrid = [];
Res.dynErr_PPA2_algoGrid = [];
Res.dynErr_PPA3_algoGrid = [];
if case_nr==2 || case_nr>3
   Res.dynErr_PPA4_algoGrid = [];
end
Res.dynErr_PFI0_algoGrid = [];
Res.dynErr_PFI1_algoGrid = [];
Res.dynErr_PFI2_algoGrid = [];
Res.dynErr_PFI3_algoGrid = [];
if case_nr==2 || case_nr>3
   Res.dynErr_PFI4_algoGrid = [];
end

%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------
for i=1:Res.nsamples
    % Euler Equ. Errors for PPA____________________________________________
    for method=0:1 
    for order_ct=0:M
        if method==0
            Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PPA.mat'));
            Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PPA.mat'));
        else
            Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
            Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'));
        end
        
        % Algorithm-implied policies
        K(method+1,order_ct+1,i) = StaticParams.kGrid*squeeze(pdf(method+1,order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
        sd(method+1,order_ct+1,i) = sqrt(StaticParams.kGrid.^2*squeeze(pdf(method+1,order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]);
        skew(method+1,order_ct+1,i) = ((StaticParams.kGrid-squeeze(K(method+1,order_ct+1,i))).^3 ...
                                      *squeeze(pdf(method+1,order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                                      ./squeeze(sd(method+1,order_ct+1,i))^3;
        kurt(method+1,order_ct+1,i) = ((StaticParams.kGrid-squeeze(K(method+1,order_ct+1,i))).^4 ...
                                      *squeeze(pdf(method+1,order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                                      ./squeeze(sd(method+1,order_ct+1,i))^4;
        gini(method+1,order_ct+1,i) = sum(sum((repmat(StaticParams.kGrid',[1,length(StaticParams.kGrid)])-repmat(StaticParams.kGrid,[length(StaticParams.kGrid),1]))...
                                      .*repmat((squeeze(pdf(method+1,order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]),[1,length(StaticParams.kGrid)])...
                                      .*repmat((squeeze(pdf(method+1,order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])',[length(StaticParams.kGrid),1]),2),1)...
                                      ./(2*squeeze(K(method+1,order_ct+1,i)));
        k_prime = interp1(StaticParams.kGrid_pol,squeeze(...
                  [getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))')])',Res.kGrid)';
        c_prime = interp1(StaticParams.kGrid_pol,squeeze(...
                  [getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))')])',Res.kGrid)';   
        % standard Euler policies
        euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                  *sum(interp1(Res.kGrid,sort(k_prime,2)',StaticParams.kGrid)'.*squeeze(pdf(method+1,order_ct+1,:,:)),2);
        euler_weights = calcWeights(squeeze(cdf(method+1,order_ct+1,:,:)),Poly,StaticParams,Res.z_ag(i),...
                            interp1(Res.kGrid,sort(k_prime,2)',StaticParams.kGrid)');
        if i<Res.nsamples
            LoM = calcAggregates(interp1(Res.kGrid,sort(k_prime,2)',StaticParams.kGrid)',...
                  squeeze(sum(repmat(squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1))),[1,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:))]).*Poly.total,1)),...
                  StaticParams,Poly);
            K_pr_LoM(method+1,order_ct+1,i) = LoM(Res.z_ag(i)+1);
        end
        euler_expec = @(grid) repmat([1-StaticParams.delta+rate(0,euler_K);...
                                      1-StaticParams.delta+rate(1,euler_K)],[1,length(grid)])...
                      .*(interp1(StaticParams.kGrid_pol,...
                               [getCApprox(Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')')...
                      .^(-StaticParams.gamma);
        euler_c_prime = min(wealth(Res.z_ag(i),squeeze(K(method+1,order_ct+1,i)),Res.kGrid),...
                                    [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime(1,:))).^(-1/StaticParams.gamma);...
                                     (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime(2,:))).^(-1/StaticParams.gamma)]);
        % standard error
        Err = (c_prime-euler_c_prime)./max(eps,euler_c_prime);
        Err_algoGrid = (c_prime(:,algo_idx)-euler_c_prime(:,algo_idx))./max(eps,euler_c_prime(:,algo_idx));
        Err = Err(logical(squeeze(idx_posProb(method+1,order_ct+1,:,:))));
        Err_algoGrid = Err_algoGrid(logical(squeeze(idx_posProb(method+1,order_ct+1,:,algo_idx))));
        switch method
            case 0
                switch order_ct
                    case 0
                    Res.stErr_PPA0 = [Res.stErr_PPA0;Err];
                    Res.stErr_PPA0_algoGrid = [Res.stErr_PPA0_algoGrid;Err_algoGrid];
                    case 1
                    Res.stErr_PPA1 = [Res.stErr_PPA1;Err];
                    Res.stErr_PPA1_algoGrid = [Res.stErr_PPA1_algoGrid;Err_algoGrid];
                    case 2
                    Res.stErr_PPA2 = [Res.stErr_PPA2;Err];
                    Res.stErr_PPA2_algoGrid = [Res.stErr_PPA2_algoGrid;Err_algoGrid];
                    case 3
                    Res.stErr_PPA3 = [Res.stErr_PPA3;Err];
                    Res.stErr_PPA3_algoGrid = [Res.stErr_PPA3_algoGrid;Err_algoGrid];
                    case 4
                    Res.stErr_PPA4 = [Res.stErr_PPA4;Err];
                    Res.stErr_PPA4_algoGrid = [Res.stErr_PPA4_algoGrid;Err_algoGrid];
                end
            case 1
                switch order_ct
                    case 0
                    Res.stErr_PFI0 = [Res.stErr_PFI0;Err];
                    Res.stErr_PFI0_algoGrid = [Res.stErr_PFI0_algoGrid;Err_algoGrid];
                    case 1
                    Res.stErr_PFI1 = [Res.stErr_PFI1;Err];
                    Res.stErr_PFI1_algoGrid = [Res.stErr_PFI1_algoGrid;Err_algoGrid];
                    case 2
                    Res.stErr_PFI2 = [Res.stErr_PFI2;Err];
                    Res.stErr_PFI2_algoGrid = [Res.stErr_PFI2_algoGrid;Err_algoGrid];
                    case 3
                    Res.stErr_PFI3 = [Res.stErr_PFI3;Err];
                    Res.stErr_PFI3_algoGrid = [Res.stErr_PFI3_algoGrid;Err_algoGrid];
                    case 4
                    Res.stErr_PFI4 = [Res.stErr_PFI4;Err];
                    Res.stErr_PFI4_algoGrid = [Res.stErr_PFI4_algoGrid;Err_algoGrid];
                end
        end
        % update distribution
        if i<Res.nsamples
            pdf(method+1,order_ct+1,:,:) = calcEoPDistr(interp1(Res.kGrid,sort(k_prime,2)',StaticParams.kGrid)',...
                           squeeze(cdf(method+1,order_ct+1,:,:)),StaticParams.kGrid);
            pdf(method+1,order_ct+1,:,:) = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                           .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                           *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(pdf,4)])...
                           .*squeeze(pdf(method+1,order_ct+1,:,:)));
            idx = max(find(squeeze(pdf(method+1,order_ct+1,1,:))>eps,1,'last'),find(squeeze(pdf(method+1,order_ct+1,2,:))>eps,1,'last'));
            idx_posProb(method+1,order_ct+1,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
            pdf = pdf.*(pdf>eps);
            pdf = pdf./repmat(sum(pdf,4),[1,1,1,size(pdf,4)]);
            cdf = cummax(min(1,cumsum(pdf,4)),4);
            cdf(method+1,order_ct+1,1,cdf(method+1,order_ct+1,1,:)==max(cdf(method+1,order_ct+1,1,:))) = 1;
            cdf(method+1,order_ct+1,2,cdf(method+1,order_ct+1,2,:)==max(cdf(method+1,order_ct+1,2,:))) = 1;
            weights(method+1,order_ct+1,1:(Poly.phi_N+1))...
            = calcWeights(squeeze(cdf(method+1,order_ct+1,:,:)),Poly,StaticParams,Res.z_ag(i+1));
        end
        
        % dynamic Euler policies and distributions
        dyn_euler_K(method+1,order_ct+1,i) = StaticParams.kGrid*squeeze(dyn_pdf(method+1,order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
        k_prime_next = [getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(dyn_weights(method+1,order_ct+1,1:(Poly.phi_N+1)))');...
                        getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(dyn_weights(method+1,order_ct+1,1:(Poly.phi_N+1)))')];  
        euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                  *sum(interp1(StaticParams.kGrid_pol,sort(k_prime_next,2)',StaticParams.kGrid)'.*squeeze(dyn_pdf(method+1,order_ct+1,:,:)),2);
        euler_weights = calcWeights(squeeze(dyn_cdf(method+1,order_ct+1,:,:)),Poly,StaticParams,Res.z_ag(i),...
                        interp1(StaticParams.kGrid_pol,sort(k_prime_next,2)',StaticParams.kGrid)');
        euler_expec = @(grid) repmat([1-StaticParams.delta+rate(0,euler_K);...
                                      1-StaticParams.delta+rate(1,euler_K)],[1,length(grid)])...
                      .*(interp1(StaticParams.kGrid_pol,...
                               [getCApprox(Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')')...
                      .^(-StaticParams.gamma);
        dyn_euler_c_prime = min(wealth(Res.z_ag(i),squeeze(dyn_euler_K(method+1,order_ct+1,i)),Res.kGrid),...
                                    [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(interp1(StaticParams.kGrid_pol,sort(k_prime_next(1,:),2)',Res.kGrid)')).^(-1/StaticParams.gamma);...
                                     (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(interp1(StaticParams.kGrid_pol,sort(k_prime_next(2,:),2)',Res.kGrid)')).^(-1/StaticParams.gamma)]);
        dyn_euler_k_prime = wealth(Res.z_ag(i),squeeze(dyn_euler_K(method+1,order_ct+1,i)),Res.kGrid)-dyn_euler_c_prime;
        % dynamic error
        Err = (c_prime-dyn_euler_c_prime)./max(eps,dyn_euler_c_prime);
        Err_algoGrid = (c_prime(:,algo_idx)-dyn_euler_c_prime(:,algo_idx))./max(eps,dyn_euler_c_prime(:,algo_idx));
        Err = Err(logical(squeeze(dyn_idx_posProb(method+1,order_ct+1,:,:))));
        Err_algoGrid = Err_algoGrid(logical(squeeze(dyn_idx_posProb(method+1,order_ct+1,:,algo_idx))));
        switch method
            case 0
                switch order_ct
                    case 0
                    Res.dynErr_PPA0 = [Res.dynErr_PPA0;Err];
                    Res.dynErr_PPA0_algoGrid = [Res.dynErr_PPA0_algoGrid;Err_algoGrid];
                    case 1
                    Res.dynErr_PPA1 = [Res.dynErr_PPA1;Err];
                    Res.dynErr_PPA1_algoGrid = [Res.dynErr_PPA1_algoGrid;Err_algoGrid];
                    case 2
                    Res.dynErr_PPA2 = [Res.dynErr_PPA2;Err];
                    Res.dynErr_PPA2_algoGrid = [Res.dynErr_PPA2_algoGrid;Err_algoGrid];
                    case 3
                    Res.dynErr_PPA3 = [Res.dynErr_PPA3;Err];
                    Res.dynErr_PPA3_algoGrid = [Res.dynErr_PPA3_algoGrid;Err_algoGrid];
                    case 4
                    Res.dynErr_PPA4 = [Res.dynErr_PPA4;Err];
                    Res.dynErr_PPA4_algoGrid = [Res.dynErr_PPA4_algoGrid;Err_algoGrid];
                end
            case 1
                switch order_ct
                    case 0
                    Res.dynErr_PFI0 = [Res.dynErr_PFI0;Err];
                    Res.dynErr_PFI0_algoGrid = [Res.dynErr_PFI0_algoGrid;Err_algoGrid];
                    case 1
                    Res.dynErr_PFI1 = [Res.dynErr_PFI1;Err];
                    Res.dynErr_PFI1_algoGrid = [Res.dynErr_PFI1_algoGrid;Err_algoGrid];
                    case 2
                    Res.dynErr_PFI2 = [Res.dynErr_PFI2;Err];
                    Res.dynErr_PFI2_algoGrid = [Res.dynErr_PFI2_algoGrid;Err_algoGrid];
                    case 3
                    Res.dynErr_PFI3 = [Res.dynErr_PFI3;Err];
                    Res.dynErr_PFI3_algoGrid = [Res.dynErr_PFI3_algoGrid;Err_algoGrid];
                    case 4
                    Res.dynErr_PFI4 = [Res.dynErr_PFI4;Err];
                    Res.dynErr_PFI4_algoGrid = [Res.dynErr_PFI4_algoGrid;Err_algoGrid];
                end
        end
        % update distribution
        if i<Res.nsamples
            dyn_pdf(method+1,order_ct+1,:,:) = calcEoPDistr(interp1(Res.kGrid,sort(dyn_euler_k_prime,2)',StaticParams.kGrid)',...
                           squeeze(dyn_cdf(method+1,order_ct+1,:,:)),StaticParams.kGrid);
            dyn_pdf(method+1,order_ct+1,:,:) = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                           .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                           *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(dyn_pdf,4)])...
                           .*squeeze(dyn_pdf(method+1,order_ct+1,:,:)));
            idx = max(find(squeeze(dyn_pdf(method+1,order_ct+1,1,:))>eps,1,'last'),find(squeeze(dyn_pdf(method+1,order_ct+1,2,:))>eps,1,'last'));
            dyn_idx_posProb(method+1,order_ct+1,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
            dyn_pdf = dyn_pdf.*(dyn_pdf>eps);
            dyn_pdf = dyn_pdf./repmat(sum(dyn_pdf,4),[1,1,1,size(dyn_pdf,4)]);
            dyn_cdf = cummax(min(1,cumsum(dyn_pdf,4)),4);
            dyn_cdf(method+1,order_ct+1,1,dyn_cdf(method+1,order_ct+1,1,:)==max(dyn_cdf(method+1,order_ct+1,1,:))) = 1;
            dyn_cdf(method+1,order_ct+1,2,dyn_cdf(method+1,order_ct+1,2,:)==max(dyn_cdf(method+1,order_ct+1,2,:))) = 1;
            dyn_weights(method+1,order_ct+1,1:(Poly.phi_N+1))...
            = calcWeights(squeeze(dyn_cdf(method+1,order_ct+1,:,:)),Poly,StaticParams,Res.z_ag(i+1));
        end
    end
    end
    
    % Euler Equ. Errors for Krusell_Smith__________________________________
    for order_ct=1:M_KS
    KS = load(strcat(folder,'/KS',num2str(order_ct),'_Sol',num2str(case_nr),'.mat'));
    % Algorithm-implied policies and distributions
    K_KS(order_ct,i) = StaticParams.kGrid*squeeze(pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
    sd_KS(order_ct,i) = sqrt(StaticParams.kGrid.^2*squeeze(pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]);
    skew_KS(order_ct,i) = ((StaticParams.kGrid-squeeze(K_KS(order_ct,i))).^3 ...
                                  *squeeze(pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                                  ./squeeze(sd_KS(order_ct,i))^3;
    kurt_KS(order_ct,i) = ((StaticParams.kGrid-squeeze(K_KS(order_ct,i))).^4 ...
                                  *squeeze(pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                                  ./squeeze(sd_KS(order_ct,i))^4;
    gini_KS(order_ct,i) = sum(sum((repmat(StaticParams.kGrid',[1,length(StaticParams.kGrid)])-repmat(StaticParams.kGrid,[length(StaticParams.kGrid),1]))...
                                  .*repmat((squeeze(pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]),[1,length(StaticParams.kGrid)])...
                                  .*repmat((squeeze(pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])',[length(StaticParams.kGrid),1]),2),1)...
                                  ./(2*squeeze(K_KS(order_ct,i)));
    k_prime_KS = [interpn(KS.k,KS.km,KS.kprime(:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))))';...
                  interpn(KS.k,KS.km,KS.kprime(:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))))'];       
    c_prime_KS = max(0,wealth(Res.z_ag(i),K_KS(order_ct,i),Res.kGrid)-k_prime_KS);
    % standard Euler policies 
    euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
              *sum(interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'.*squeeze(pdf_KS(order_ct,:,:)),2);
    K_pr_LoM_KS(order_ct,i) = exp(KS.B1(2*Res.z_ag(i)+1)+KS.B1(2*Res.z_ag(i)+2)*log(K_KS(order_ct,i)));
    euler_expec = @(grid) repmat([1-StaticParams.delta+rate(0,euler_K);...
                                  1-StaticParams.delta+rate(1,euler_K)],[1,length(grid)])...
                  .*max(1e-15,([wealth(0,euler_K,grid);wealth(1,euler_K,grid)]...
                        -interp1(StaticParams.kGrid,...
                          [interpn(KS.k,KS.km,KS.kprime(:,:,1,1),StaticParams.kGrid,euler_K)';...
                           interpn(KS.k,KS.km,KS.kprime(:,:,1,2),StaticParams.kGrid,euler_K)';...
                           interpn(KS.k,KS.km,KS.kprime(:,:,2,1),StaticParams.kGrid,euler_K)';...
                           interpn(KS.k,KS.km,KS.kprime(:,:,2,2),StaticParams.kGrid,euler_K)']',...
                         grid,'linear','extrap')'))...
                  .^(-StaticParams.gamma);
    euler_c_prime_KS = min(wealth(Res.z_ag(i),K_KS(order_ct,i),Res.kGrid),...
                              [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime_KS(1,:))).^(-1/StaticParams.gamma);...
                               (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime_KS(2,:))).^(-1/StaticParams.gamma)]);
    % standard error
    Err = (c_prime_KS-euler_c_prime_KS)./max(eps,euler_c_prime_KS);
    Err_algoGrid = (c_prime_KS(:,algo_idx)-euler_c_prime_KS(:,algo_idx))./max(eps,euler_c_prime_KS(:,algo_idx));
    Err = Err(logical(squeeze(idx_posProb_KS(order_ct,:,:))));
    Err_algoGrid = Err_algoGrid(logical(squeeze(idx_posProb_KS(order_ct,:,algo_idx))));
    switch order_ct
        case 1
        Res.stErr_KS1 = [Res.stErr_KS1;Err];
        Res.stErr_KS1_algoGrid = [Res.stErr_KS1_algoGrid;Err_algoGrid];
        case 2
        Res.stErr_KS2 = [Res.stErr_KS2;Err];
        Res.stErr_KS2_algoGrid = [Res.stErr_KS2_algoGrid;Err_algoGrid];
        case 3
        Res.stErr_KS3 = [Res.stErr_KS3;Err];
        Res.stErr_KS3_algoGrid = [Res.stErr_KS3_algoGrid;Err_algoGrid];
        case 4
        Res.stErr_KS4 = [Res.stErr_KS4;Err];
        Res.stErr_KS4_algoGrid = [Res.stErr_KS4_algoGrid;Err_algoGrid];
    end
    % update distribution
    if i<Res.nsamples
        pdf_KS(order_ct,:,:) = calcEoPDistr(interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)',...
                 squeeze(cdf_KS(order_ct,:,:)),StaticParams.kGrid);        
        pdf_KS(order_ct,:,:) = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                   .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                   *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(pdf_KS,3)])...
                   .*squeeze(pdf_KS(order_ct,:,:)));
        idx = max(find(squeeze(pdf_KS(order_ct,1,:))>eps,1,'last'),find(squeeze(pdf_KS(order_ct,2,:))>eps,1,'last'));
        idx_posProb_KS(:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
        pdf_KS = pdf_KS.*(pdf_KS>eps);
        pdf_KS = pdf_KS./repmat(sum(pdf_KS,3),[1,1,size(pdf_KS,3)]);   
        cdf_KS = cummax(min(1,cumsum(pdf_KS,3)),3);
        cdf_KS(order_ct,1,cdf_KS(order_ct,1,:)==max(cdf_KS(order_ct,1,:))) = 1;
        cdf_KS(order_ct,2,cdf_KS(order_ct,2,:)==max(cdf_KS(order_ct,2,:))) = 1; 
    end
    
    % dynamic Euler policies and distributions
    dyn_euler_K_KS(order_ct,i) = StaticParams.kGrid*squeeze(dyn_pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
    k_prime_next = [interpn(KS.k,KS.km,KS.kprime(:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))))';...
                    interpn(KS.k,KS.km,KS.kprime(:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))))'];       
    euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
              *sum(interp1(Res.kGrid,sort(k_prime_next,2)',StaticParams.kGrid)'.*squeeze(dyn_pdf_KS(order_ct,:,:)),2);
    euler_expec = @(grid) repmat([1-StaticParams.delta+rate(0,euler_K);...
                                  1-StaticParams.delta+rate(1,euler_K)],[1,length(grid)])...
                  .*max(1e-15,([wealth(0,euler_K,grid);wealth(1,euler_K,grid)]...
                        -interp1(StaticParams.kGrid,...
                          [interpn(KS.k,KS.km,KS.kprime(:,:,1,1),StaticParams.kGrid,euler_K)';...
                           interpn(KS.k,KS.km,KS.kprime(:,:,1,2),StaticParams.kGrid,euler_K)';...
                           interpn(KS.k,KS.km,KS.kprime(:,:,2,1),StaticParams.kGrid,euler_K)';...
                           interpn(KS.k,KS.km,KS.kprime(:,:,2,2),StaticParams.kGrid,euler_K)']',...
                         grid','linear','extrap')'))...
                  .^(-StaticParams.gamma);
    dyn_euler_c_prime_KS = min(wealth(Res.z_ag(i),dyn_euler_K_KS(order_ct,i),Res.kGrid),...
                              [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime_next(1,:))).^(-1/StaticParams.gamma);...
                               (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime_next(2,:))).^(-1/StaticParams.gamma)]);
    dyn_euler_k_prime_KS = wealth(Res.z_ag(i),dyn_euler_K_KS(order_ct,i),Res.kGrid)-dyn_euler_c_prime_KS;
    % dynamic error
    Err = (c_prime_KS-dyn_euler_c_prime_KS)./max(eps,dyn_euler_c_prime_KS);
    Err_algoGrid = (c_prime_KS(:,algo_idx)-dyn_euler_c_prime_KS(:,algo_idx))./max(eps,dyn_euler_c_prime_KS(:,algo_idx));
    Err = Err(logical(squeeze(dyn_idx_posProb_KS(order_ct,:,:))));
    Err_algoGrid = Err_algoGrid(logical(squeeze(dyn_idx_posProb_KS(order_ct,:,algo_idx))));
    switch order_ct
        case 1
        Res.dynErr_KS1 = [Res.dynErr_KS1;Err];
        Res.dynErr_KS1_algoGrid = [Res.dynErr_KS1_algoGrid;Err_algoGrid];
        case 2
        Res.dynErr_KS2 = [Res.dynErr_KS2;Err];
        Res.dynErr_KS2_algoGrid = [Res.dynErr_KS2_algoGrid;Err_algoGrid];
        case 3
        Res.dynErr_KS3 = [Res.dynErr_KS3;Err];
        Res.dynErr_KS3_algoGrid = [Res.dynErr_KS3_algoGrid;Err_algoGrid];
        case 4
        Res.dynErr_KS4 = [Res.dynErr_KS4;Err];
        Res.dynErr_KS4_algoGrid = [Res.dynErr_KS4_algoGrid;Err_algoGrid];
    end
    % update distribution
    if i<Res.nsamples
        dyn_pdf_KS(order_ct,:,:) = calcEoPDistr(interp1(Res.kGrid,sort(dyn_euler_k_prime_KS,2)',StaticParams.kGrid)',...
                     squeeze(dyn_cdf_KS(order_ct,:,:)),StaticParams.kGrid);
        dyn_pdf_KS(order_ct,:,:) = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                   .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                   *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(dyn_pdf_KS,3)])...
                   .*squeeze(dyn_pdf_KS(order_ct,:,:)));  
        idx = max(find(squeeze(dyn_pdf_KS(order_ct,1,:))>eps,1,'last'),find(squeeze(dyn_pdf_KS(order_ct,2,:))>eps,1,'last'));
        dyn_idx_posProb_KS(:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
        dyn_pdf_KS = dyn_pdf_KS.*(dyn_pdf_KS>eps);
        dyn_pdf_KS = dyn_pdf_KS./repmat(sum(dyn_pdf_KS,3),[1,1,size(dyn_pdf_KS,3)]);
        dyn_cdf_KS = cummax(min(1,cumsum(dyn_pdf_KS,3)),3);
        dyn_cdf_KS(order_ct,1,dyn_cdf_KS(order_ct,1,:)==max(dyn_cdf_KS(order_ct,1,:))) = 1;
        dyn_cdf_KS(order_ct,2,dyn_cdf_KS(order_ct,2,:)==max(dyn_cdf_KS(order_ct,2,:))) = 1;  
    end
    end
    
    % Euler Equ. Errors for Reiter_________________________________________
    % Algorithm-implied policies and distributions
    K_R(i) = StaticParams.kGrid*pdf_R'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
    sd_R(i) = sqrt(StaticParams.kGrid.^2*pdf_R'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]);
    skew_R(i) = ((StaticParams.kGrid-K_R(i)).^3 ...
                  *pdf_R'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                  ./sd_R(i)^3;
    kurt_R(i) = ((StaticParams.kGrid-K_R(i)).^4 ...
                  *pdf_R'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                  ./sd_R(i)^4;
    gini_R(i) = sum(sum((repmat(StaticParams.kGrid',[1,length(StaticParams.kGrid)])-repmat(StaticParams.kGrid,[length(StaticParams.kGrid),1]))...
                  .*repmat((pdf_R'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]),[1,length(StaticParams.kGrid)])...
                  .*repmat((pdf_R'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])',[length(StaticParams.kGrid),1]),2),1)...
                  ./(2*K_R(i));
    k_prime_R = interp1(StaticParams.kGrid_pol,squeeze(...
                           [getCApprox(R.Params2.GridMat,squeeze(R.k_prime(:,2*Res.z_ag(i)+1,:)),R.idx,K_R(i));...
                            getCApprox(R.Params2.GridMat,squeeze(R.k_prime(:,2*Res.z_ag(i)+2,:)),R.idx,K_R(i))]...       
                            )',Res.kGrid)';
	c_prime_R = interp1(StaticParams.kGrid_pol,squeeze(...
                           [getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,2*Res.z_ag(i)+1,:)),R.idx,K_R(i));...
                            getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,2*Res.z_ag(i)+2,:)),R.idx,K_R(i))]...       
                            )',Res.kGrid)';   
    % standard Euler policies
    euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
              *sum(interp1(Res.kGrid,sort(k_prime_R,2)',StaticParams.kGrid)'.*pdf_R,2);
    euler_expec = @(grid) repmat([1-StaticParams.delta+rate(0,euler_K);...
                                  1-StaticParams.delta+rate(1,euler_K)],[1,length(grid)])...
                  .*(interp1(StaticParams.kGrid_pol,...
                          [getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,1,:)),R.idx,euler_K);...
                           getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,2,:)),R.idx,euler_K);...
                           getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,3,:)),R.idx,euler_K);...
                           getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,4,:)),R.idx,euler_K)]',...
                         grid','linear','extrap')')...
                  .^(-StaticParams.gamma);
    euler_c_prime_R = min(wealth(Res.z_ag(i),K_R(i),Res.kGrid),...
                              [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime_R(1,:))).^(-1/StaticParams.gamma);...
                               (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime_R(2,:))).^(-1/StaticParams.gamma)]);
    % standard error
    Err = (c_prime_R-euler_c_prime_R)./max(eps,euler_c_prime_R);
    Err_algoGrid = (c_prime_R(:,algo_idx)-euler_c_prime_R(:,algo_idx))./max(eps,euler_c_prime_R(:,algo_idx));
    Err = Err(logical(idx_posProb_R));
    Err_algoGrid = Err_algoGrid(logical(idx_posProb_R(:,algo_idx)));
    Res.stErr_R = [Res.stErr_R;Err];
    Res.stErr_R_algoGrid = [Res.stErr_R_algoGrid;Err_algoGrid];
    % update distribution
    if i<Res.nsamples
        pdf_R = calcEoPDistr(interp1(Res.kGrid,sort(k_prime_R,2)',StaticParams.kGrid)',...
                 cdf_R,StaticParams.kGrid);
        pdf_R = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                   .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                   *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(pdf_R,2)])...
                   .*pdf_R);
        idx = max(find(squeeze(pdf_R(1,:))>eps,1,'last'),find(squeeze(pdf_R(2,:))>eps,1,'last'));
        idx_posProb_R(:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
        pdf_R = pdf_R.*(pdf_R>eps);
        pdf_R = pdf_R./repmat(sum(pdf_R,2),[1,size(pdf_R,2)]);
        cdf_R = cummax(min(1,cumsum(pdf_R,2)),2);
        cdf_R(1,cdf_R(1,:)==max(cdf_R(1,:))) = 1;
        cdf_R(2,cdf_R(2,:)==max(cdf_R(2,:))) = 1;
    end
    
    % dynamic Euler policies and distributions
    dyn_euler_K_R(i) = StaticParams.kGrid*dyn_pdf_R'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
    k_prime_next = [getCApprox(R.Params2.GridMat,squeeze(R.k_prime(:,2*Res.z_ag(i)+1,:)),R.idx,dyn_euler_K_R(i));...
                    getCApprox(R.Params2.GridMat,squeeze(R.k_prime(:,2*Res.z_ag(i)+2,:)),R.idx,dyn_euler_K_R(i))];       
    euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
              *sum(interp1(StaticParams.kGrid_pol,sort(k_prime_next,2)',StaticParams.kGrid)'.*dyn_pdf_R,2);
    euler_expec = @(grid) repmat([1-StaticParams.delta+rate(0,euler_K);...
                                  1-StaticParams.delta+rate(1,euler_K)],[1,length(grid)])...
                  .*(interp1(StaticParams.kGrid_pol,...
                          [getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,1,:)),R.idx,euler_K);...
                           getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,2,:)),R.idx,euler_K);...
                           getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,3,:)),R.idx,euler_K);...
                           getCApprox(R.Params2.GridMat,squeeze(R.c_prime(:,4,:)),R.idx,euler_K)]',...
                         grid,'linear','extrap')')...
                  .^(-StaticParams.gamma);
    dyn_euler_c_prime_R = min(wealth(Res.z_ag(i),dyn_euler_K_R(i),Res.kGrid),...
                              [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(interp1(StaticParams.kGrid_pol,sort(k_prime_next(1,:),2)',Res.kGrid)')).^(-1/StaticParams.gamma);...
                               (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(interp1(StaticParams.kGrid_pol,sort(k_prime_next(2,:),2)',Res.kGrid)')).^(-1/StaticParams.gamma)]);
    dyn_euler_k_prime_R = wealth(Res.z_ag(i),dyn_euler_K_R(i),Res.kGrid)-dyn_euler_c_prime_R;
    % dynamic error
    Err = (c_prime_R-dyn_euler_c_prime_R)./max(eps,dyn_euler_c_prime_R);
    Err_algoGrid = (c_prime_R(:,algo_idx)-dyn_euler_c_prime_R(:,algo_idx))./max(eps,dyn_euler_c_prime_R(:,algo_idx));
    Err = Err(logical(dyn_idx_posProb_R));
    Err_algoGrid = Err_algoGrid(logical(dyn_idx_posProb_R(:,algo_idx)));
    Res.dynErr_R = [Res.dynErr_R;Err];
    Res.dynErr_R_algoGrid = [Res.dynErr_R_algoGrid;Err_algoGrid];
    % update distribution
    if i<Res.nsamples
        dyn_pdf_R = calcEoPDistr(interp1(Res.kGrid,sort(dyn_euler_k_prime_R,2)',StaticParams.kGrid)',...
                     dyn_cdf_R,StaticParams.kGrid);
        dyn_pdf_R = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                   .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                   *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(dyn_pdf_R,2)])...
                   .*dyn_pdf_R);
        idx = max(find(squeeze(dyn_pdf_R(1,:))>eps,1,'last'),find(squeeze(dyn_pdf_R(2,:))>eps,1,'last'));
        dyn_idx_posProb_R(:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
        dyn_pdf_R = dyn_pdf_R.*(dyn_pdf_R>eps);
        dyn_pdf_R = dyn_pdf_R./repmat(sum(dyn_pdf_R,2),[1,size(dyn_pdf_R,2)]);
        dyn_cdf_R = cummax(min(1,cumsum(dyn_pdf_R,2)),2);
        dyn_cdf_R(1,dyn_cdf_R(1,:)==max(dyn_cdf_R(1,:))) = 1;
        dyn_cdf_R(2,dyn_cdf_R(2,:)==max(dyn_cdf_R(2,:))) = 1;
    end
    
    % Euler Equ. Errors for Den Haan-Rendahl_______________________________
    % Algorithm-implied policies and distributions
    K_DR(i) = StaticParams.kGrid*pdf_DR'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
    sd_DR(i) = sqrt(StaticParams.kGrid.^2*pdf_DR'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]);
    skew_DR(i) = ((StaticParams.kGrid-K_DR(i)).^3 ...
                  *pdf_DR'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                  ./sd_DR(i)^3;
    kurt_DR(i) = ((StaticParams.kGrid-K_DR(i)).^4 ...
                  *pdf_DR'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                  ./sd_DR(i)^4;
    gini_DR(i) = sum(sum((repmat(StaticParams.kGrid',[1,length(StaticParams.kGrid)])-repmat(StaticParams.kGrid,[length(StaticParams.kGrid),1]))...
                  .*repmat((pdf_DR'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]),[1,length(StaticParams.kGrid)])...
                  .*repmat((pdf_DR'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])',[length(StaticParams.kGrid),1]),2),1)...
                  ./(2*K_DR(i));
    k_prime_DR = [interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,2,2-Res.z_ag(i)),Res.kGrid,max(min(StaticParams.kGrid*pdf_DR(2,:)',DR.Kemax),DR.Kemin),max(min(StaticParams.kGrid*pdf_DR(1,:)',DR.Kumax),DR.Kumin))';...
                             interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,1,2-Res.z_ag(i)),Res.kGrid,max(min(StaticParams.kGrid*pdf_DR(2,:)',DR.Kemax),DR.Kemin),max(min(StaticParams.kGrid*pdf_DR(1,:)',DR.Kumax),DR.Kumin))'];       
    c_prime_DR = max(0,wealth(Res.z_ag(i),K_DR(i),Res.kGrid)-k_prime_DR);
    % standard Euler policies 
    euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
              *sum(interp1(Res.kGrid,sort(k_prime_DR,2)',StaticParams.kGrid)'.*pdf_DR,2);
    euler_K_cond = sum(interp1(Res.kGrid,sort(k_prime_DR,2)',StaticParams.kGrid)'.*pdf_DR,2);
    euler_K_cond = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),:))',[1,2])...
                   .*StaticParams.P(2*Res.z_ag(i)+(1:2),:)')...
                   *(StaticParams.p(2*Res.z_ag(i)+(1:2)).*euler_K_cond);     
    euler_LoM_empl = @(z_ag_pr_idx) interpn(DR.Ke,DR.Ku,DR.KpeN(:,:,Res.z_ag(i)+1,z_ag_pr_idx),...
                     max(min(StaticParams.kGrid*pdf_DR(2,:)',DR.Kemax),DR.Kemin),max(min(StaticParams.kGrid*pdf_DR(1,:)',DR.Kumax),DR.Kumin))';
    euler_LoM_unempl =  @(z_ag_pr_idx) interpn(DR.Ke,DR.Ku,DR.KpuN(:,:,Res.z_ag(i)+1,z_ag_pr_idx),...
                     max(min(StaticParams.kGrid*pdf_DR(2,:)',DR.Kemax),DR.Kemin),max(min(StaticParams.kGrid*pdf_DR(1,:)',DR.Kumax),DR.Kumin))';
    K_pr_LoM_DR(i) = StaticParams.p(2*Res.z_ag(i)+(1:2))'*[euler_LoM_unempl(Res.z_ag(i)+1);euler_LoM_empl(Res.z_ag(i)+1)];
    euler_expec = @(grid) repmat([1-StaticParams.delta+rate(0,euler_K);...
                                  1-StaticParams.delta+rate(1,euler_K)],[1,length(grid)])...
                  .*max(1e-15,([wealth(0,euler_K,grid);wealth(1,euler_K,grid)]...
                        -interp1(StaticParams.kGrid,...
                          [interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,2,2),StaticParams.kGrid,euler_K_cond(2),euler_K_cond(1))';...
                           interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,1,2),StaticParams.kGrid,euler_K_cond(2),euler_K_cond(1))';...
                           interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,2,1),StaticParams.kGrid,euler_K_cond(4),euler_K_cond(3))';...
                           interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,1,1),StaticParams.kGrid,euler_K_cond(4),euler_K_cond(3))']',...
                         grid','linear','extrap')'))...
                  .^(-StaticParams.gamma);
    euler_c_prime_DR = min(wealth(Res.z_ag(i),K_DR(i),Res.kGrid),...
                              [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime_DR(1,:))).^(-1/StaticParams.gamma);...
                               (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime_DR(2,:))).^(-1/StaticParams.gamma)]);
    % standard error
    Err = (c_prime_DR-euler_c_prime_DR)./max(eps,euler_c_prime_DR);
    Err_algoGrid = (c_prime_DR(:,algo_idx)-euler_c_prime_DR(:,algo_idx))./max(eps,euler_c_prime_DR(:,algo_idx));
    Err = Err(logical(idx_posProb_DR));
    Err_algoGrid = Err_algoGrid(logical(idx_posProb_DR(:,algo_idx)));
    Res.stErr_DR = [Res.stErr_DR;Err];
    Res.stErr_DR_algoGrid = [Res.stErr_DR_algoGrid;Err_algoGrid];
    % update distribution
    if i<Res.nsamples
        pdf_DR = calcEoPDistr(interp1(Res.kGrid,sort(k_prime_DR,2)',StaticParams.kGrid)',...
                 cdf_DR,StaticParams.kGrid);
        pdf_DR = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                   .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                   *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(pdf_DR,2)])...
                   .*pdf_DR);
        idx = max(find(squeeze(pdf_DR(1,:))>eps,1,'last'),find(squeeze(pdf_DR(2,:))>eps,1,'last'));
        idx_posProb_DR(:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
        pdf_DR = pdf_DR.*(pdf_DR>eps);
        pdf_DR = pdf_DR./repmat(sum(pdf_DR,2),[1,size(pdf_DR,2)]);
        cdf_DR = cummax(min(1,cumsum(pdf_DR,2)),2);
        cdf_DR(1,cdf_DR(1,:)==max(cdf_DR(1,:))) = 1;
        cdf_DR(2,cdf_DR(2,:)==max(cdf_DR(2,:))) = 1;
    end
    
    % dynamic Euler policies and distributions
    dyn_euler_K_DR(i) = StaticParams.kGrid*dyn_pdf_DR'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
    k_prime_next = [interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,2,2-Res.z_ag(i)),Res.kGrid,...
                    max(min(StaticParams.kGrid*dyn_pdf_DR(2,:)',DR.Kemax),DR.Kemin),max(min(StaticParams.kGrid*dyn_pdf_DR(1,:)',DR.Kumax),DR.Kumin))';...
                    interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,1,2-Res.z_ag(i)),Res.kGrid,...
                    max(min(StaticParams.kGrid*dyn_pdf_DR(2,:)',DR.Kemax),DR.Kemin),max(min(StaticParams.kGrid*dyn_pdf_DR(1,:)',DR.Kumax),DR.Kumin))'];             
    euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
              *sum(interp1(Res.kGrid,sort(k_prime_next,2)',StaticParams.kGrid)'.*dyn_pdf_DR,2);
    euler_K_cond = sum(interp1(Res.kGrid,sort(k_prime_next,2)',StaticParams.kGrid)'.*dyn_pdf_DR,2);
    euler_K_cond = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),:))',[1,2])...
                   .*StaticParams.P(2*Res.z_ag(i)+(1:2),:)')...
                   *(StaticParams.p(2*Res.z_ag(i)+(1:2)).*euler_K_cond);        
    euler_expec = @(grid) repmat([1-StaticParams.delta+rate(0,euler_K);...
                                  1-StaticParams.delta+rate(1,euler_K)],[1,length(grid)])...
                  .*max(1e-15,([wealth(0,euler_K,grid);wealth(1,euler_K,grid)]...
                        -interp1(StaticParams.kGrid,...
                          [interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,2,2),StaticParams.kGrid,euler_K_cond(2),euler_K_cond(1))';...
                           interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,1,2),StaticParams.kGrid,euler_K_cond(2),euler_K_cond(1))';...
                           interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,2,1),StaticParams.kGrid,euler_K_cond(4),euler_K_cond(3))';...
                           interpn(DR.kp,DR.Ke,DR.Ku,DR.kpmat(:,:,:,1,1),StaticParams.kGrid,euler_K_cond(4),euler_K_cond(3))']',...
                         grid','linear','extrap')'))...
                  .^(-StaticParams.gamma);
    dyn_euler_c_prime_DR = min(wealth(Res.z_ag(i),dyn_euler_K_DR(i),Res.kGrid),...
                              [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime_next(1,:))).^(-1/StaticParams.gamma);...
                               (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime_next(2,:))).^(-1/StaticParams.gamma)]);
    dyn_euler_k_prime_DR = wealth(Res.z_ag(i),dyn_euler_K_DR(i),Res.kGrid)-dyn_euler_c_prime_DR;
    % dynamic error
    Err = (c_prime_DR-dyn_euler_c_prime_DR)./max(eps,dyn_euler_c_prime_DR);
    Err_algoGrid = (c_prime_DR(:,algo_idx)-dyn_euler_c_prime_DR(:,algo_idx))./max(eps,dyn_euler_c_prime_DR(:,algo_idx));
    Err = Err(logical(dyn_idx_posProb_DR));
    Err_algoGrid = Err_algoGrid(logical(dyn_idx_posProb_DR(:,algo_idx)));
    Res.dynErr_DR = [Res.dynErr_DR;Err];
    Res.dynErr_DR_algoGrid = [Res.dynErr_DR_algoGrid;Err_algoGrid];
    % update distribution
    if i<Res.nsamples
        dyn_pdf_DR = calcEoPDistr(interp1(Res.kGrid,sort(dyn_euler_k_prime_DR,2)',StaticParams.kGrid)',...
                     dyn_cdf_DR,StaticParams.kGrid);
        dyn_pdf_DR = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                   .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                   *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(dyn_pdf_DR,2)])...
                   .*dyn_pdf_DR);
        idx = max(find(squeeze(dyn_pdf_DR(1,:))>eps,1,'last'),find(squeeze(dyn_pdf_DR(2,:))>eps,1,'last'));
        dyn_idx_posProb_DR(:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
        dyn_pdf_DR = dyn_pdf_DR.*(dyn_pdf_DR>eps);
        dyn_pdf_DR = dyn_pdf_DR./repmat(sum(dyn_pdf_DR,2),[1,size(dyn_pdf_DR,2)]);
        dyn_cdf_DR = cummax(min(1,cumsum(dyn_pdf_DR,2)),2);
        dyn_cdf_DR(1,dyn_cdf_DR(1,:)==max(dyn_cdf_DR(1,:))) = 1;
        dyn_cdf_DR(2,dyn_cdf_DR(2,:)==max(dyn_cdf_DR(2,:))) = 1;
    
        if Res.z_ag(i)==0
            fill(x,y,'r','EdgeColor','none')
        else
            fill(x,y,'g','EdgeColor','none')
        end
        set(gca,'nextplot','add');
        plot(StaticParams.kGrid,2*StaticParams.p(2*Res.z_ag(i+1)+(1:2))'*squeeze(pdf_KS(1,:,:)),'m','LineWidth',2);
        set(gca,'nextplot','add');
        plot(StaticParams.kGrid,2*StaticParams.p(2*Res.z_ag(i+1)+(1:2))'*pdf_R,'c','LineWidth',2);
        set(gca,'nextplot','add');
        plot(StaticParams.kGrid,2*StaticParams.p(2*Res.z_ag(i+1)+(1:2))'*squeeze(pdf(1,M,:,:)),'b','LineWidth',2);
        set(gca,'nextplot','replacechildren');
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
close(v)

Res.K_KS = K_KS;
Res.sd_KS = sd_KS;
Res.skew_KS = skew_KS;
Res.kurt_KS = kurt_KS;
Res.gini_KS = gini_KS;
Res.K_pr_LoM_KS = K_pr_LoM_KS;
Res.dyn_euler_K_KS = dyn_euler_K_KS;
Res.K_R = K_R;
Res.sd_R = sd_R;
Res.skew_R = skew_R;
Res.kurt_R = kurt_R;
Res.gini_R = gini_R;
Res.dyn_euler_K_R = dyn_euler_K_R;
Res.K_DR = K_DR;
Res.sd_DR = sd_DR;
Res.skew_DR = skew_DR;
Res.kurt_DR = kurt_DR;
Res.gini_DR = gini_DR;
Res.K_pr_LoM_DR = K_pr_LoM_DR;
Res.dyn_euler_K_DR = dyn_euler_K_DR;
Res.K = K;
Res.sd = sd;
Res.skew = skew;
Res.kurt = kurt;
Res.gini = gini;
Res.K_pr_LoM = K_pr_LoM;
Res.dyn_euler_K = dyn_euler_K;

%--------------------------------------------------------------------------
% Compute the Euler equation error statistics
%--------------------------------------------------------------------------

Res.median_stErr_KS = median(abs(Res.stErr_KS1));
if case_nr<3
    Res.median_stErr_KS = [Res.median_stErr_KS;median(abs(Res.stErr_KS2));median(abs(Res.stErr_KS3))];
end
if case_nr == 1
    Res.median_stErr_KS = [Res.median_stErr_KS;median(abs(Res.stErr_KS4))];
end
Res.median_stErr_R = median(abs(Res.stErr_R));
Res.median_stErr_DR = median(abs(Res.stErr_DR));
Res.median_stErr_PPA = [median(abs(Res.stErr_PPA0));median(abs(Res.stErr_PPA1));median(abs(Res.stErr_PPA2));median(abs(Res.stErr_PPA3))];
Res.median_stErr_PFI = [median(abs(Res.stErr_PFI0));median(abs(Res.stErr_PFI1));median(abs(Res.stErr_PFI2));median(abs(Res.stErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.median_stErr_PPA = [Res.median_stErr_PPA;median(abs(Res.stErr_PPA4))];
    Res.median_stErr_PFI = [Res.median_stErr_PFI;median(abs(Res.stErr_PFI4))];
end
Res.mean_stErr_KS = mean(abs(Res.stErr_KS1));
if case_nr<3
    Res.mean_stErr_KS = [Res.mean_stErr_KS;mean(abs(Res.stErr_KS2));mean(abs(Res.stErr_KS3))];
end
if case_nr == 1
    Res.mean_stErr_KS = [Res.mean_stErr_KS;mean(abs(Res.stErr_KS4))];
end
Res.mean_stErr_R = mean(abs(Res.stErr_R));
Res.mean_stErr_DR = mean(abs(Res.stErr_DR));
Res.mean_stErr_PPA = [mean(abs(Res.stErr_PPA0));mean(abs(Res.stErr_PPA1));mean(abs(Res.stErr_PPA2));mean(abs(Res.stErr_PPA3))];
Res.mean_stErr_PFI = [mean(abs(Res.stErr_PFI0));mean(abs(Res.stErr_PFI1));mean(abs(Res.stErr_PFI2));mean(abs(Res.stErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.mean_stErr_PPA = [Res.mean_stErr_PPA;mean(abs(Res.stErr_PPA4))];
    Res.mean_stErr_PFI = [Res.mean_stErr_PFI;mean(abs(Res.stErr_PFI4))];
end
Res.min_stErr_KS = min(abs(Res.stErr_KS1));
if case_nr<3
    Res.min_stErr_KS = [Res.min_stErr_KS;min(abs(Res.stErr_KS2));min(abs(Res.stErr_KS3))];
end
if case_nr == 1
    Res.min_stErr_KS = [Res.min_stErr_KS;min(abs(Res.stErr_KS4))];
end
Res.min_stErr_R = min(abs(Res.stErr_R));
Res.min_stErr_DR = min(abs(Res.stErr_DR));
Res.min_stErr_PPA = [min(abs(Res.stErr_PPA0));min(abs(Res.stErr_PPA1));min(abs(Res.stErr_PPA2));min(abs(Res.stErr_PPA3))];
Res.min_stErr_PFI = [min(abs(Res.stErr_PFI0));min(abs(Res.stErr_PFI1));min(abs(Res.stErr_PFI2));min(abs(Res.stErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.min_stErr_PPA = [Res.min_stErr_PPA;min(abs(Res.stErr_PPA4))];
    Res.min_stErr_PFI = [Res.min_stErr_PFI;min(abs(Res.stErr_PFI4))];
end
Res.max_stErr_KS = max(abs(Res.stErr_KS1));
if case_nr<3
    Res.max_stErr_KS = [Res.max_stErr_KS;max(abs(Res.stErr_KS2));max(abs(Res.stErr_KS3))];
end
if case_nr == 1
    Res.max_stErr_KS = [Res.max_stErr_KS;max(abs(Res.stErr_KS4))];
end
Res.max_stErr_R = max(abs(Res.stErr_R));
Res.max_stErr_DR = max(abs(Res.stErr_DR));
Res.max_stErr_PPA = [max(abs(Res.stErr_PPA0));max(abs(Res.stErr_PPA1));max(abs(Res.stErr_PPA2));max(abs(Res.stErr_PPA3))];
Res.max_stErr_PFI = [max(abs(Res.stErr_PFI0));max(abs(Res.stErr_PFI1));max(abs(Res.stErr_PFI2));max(abs(Res.stErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.max_stErr_PPA = [Res.max_stErr_PPA;max(abs(Res.stErr_PPA4))];
    Res.max_stErr_PFI = [Res.max_stErr_PFI;max(abs(Res.stErr_PFI4))];
end

Res.median_dynErr_KS = median(abs(Res.dynErr_KS1));
if case_nr<3
    Res.median_dynErr_KS = [Res.median_dynErr_KS;median(abs(Res.dynErr_KS2));median(abs(Res.dynErr_KS3))];
end
if case_nr == 1
    Res.median_dynErr_KS = [Res.median_dynErr_KS;median(abs(Res.dynErr_KS4))];
end
Res.median_dynErr_R = median(abs(Res.dynErr_R));
Res.median_dynErr_DR = median(abs(Res.dynErr_DR));
Res.median_dynErr_PPA = [median(abs(Res.dynErr_PPA0));median(abs(Res.dynErr_PPA1));median(abs(Res.dynErr_PPA2));median(abs(Res.dynErr_PPA3))];
Res.median_dynErr_PFI = [median(abs(Res.dynErr_PFI0));median(abs(Res.dynErr_PFI1));median(abs(Res.dynErr_PFI2));median(abs(Res.dynErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.median_dynErr_PPA = [Res.median_dynErr_PPA;median(abs(Res.dynErr_PPA4))];
    Res.median_dynErr_PFI = [Res.median_dynErr_PFI;median(abs(Res.dynErr_PFI4))];
end
Res.mean_dynErr_KS = mean(abs(Res.dynErr_KS1));
if case_nr<3
    Res.mean_dynErr_KS = [Res.mean_dynErr_KS;mean(abs(Res.dynErr_KS2));mean(abs(Res.dynErr_KS3))];
end
if case_nr == 1
    Res.mean_dynErr_KS = [Res.mean_dynErr_KS;mean(abs(Res.dynErr_KS4))];
end
Res.mean_dynErr_R = mean(abs(Res.dynErr_R));
Res.mean_dynErr_DR = mean(abs(Res.dynErr_DR));
Res.mean_dynErr_PPA = [mean(abs(Res.dynErr_PPA0));mean(abs(Res.dynErr_PPA1));mean(abs(Res.dynErr_PPA2));mean(abs(Res.dynErr_PPA3))];
Res.mean_dynErr_PFI = [mean(abs(Res.dynErr_PFI0));mean(abs(Res.dynErr_PFI1));mean(abs(Res.dynErr_PFI2));mean(abs(Res.dynErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.mean_dynErr_PPA = [Res.mean_dynErr_PPA;mean(abs(Res.dynErr_PPA4))];
    Res.mean_dynErr_PFI = [Res.mean_dynErr_PFI;mean(abs(Res.dynErr_PFI4))];
end
Res.min_dynErr_KS = min(abs(Res.dynErr_KS1));
if case_nr<3
    Res.min_dynErr_KS = [Res.min_dynErr_KS;min(abs(Res.dynErr_KS2));min(abs(Res.dynErr_KS3))];
end
if case_nr == 1
    Res.min_dynErr_KS = [Res.min_dynErr_KS;min(abs(Res.dynErr_KS4))];
end
Res.min_dynErr_R = min(abs(Res.dynErr_R));
Res.min_dynErr_DR = min(abs(Res.dynErr_DR));
Res.min_dynErr_PPA = [min(abs(Res.dynErr_PPA0));min(abs(Res.dynErr_PPA1));min(abs(Res.dynErr_PPA2));min(abs(Res.dynErr_PPA3))];
Res.min_dynErr_PFI = [min(abs(Res.dynErr_PFI0));min(abs(Res.dynErr_PFI1));min(abs(Res.dynErr_PFI2));min(abs(Res.dynErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.min_dynErr_PPA = [Res.min_dynErr_PPA;min(abs(Res.dynErr_PPA4))];
    Res.min_dynErr_PFI = [Res.min_dynErr_PFI;min(abs(Res.dynErr_PFI4))];
end
Res.max_dynErr_KS = max(abs(Res.dynErr_KS1));
if case_nr<3
    Res.max_dynErr_KS = [Res.max_dynErr_KS;max(abs(Res.dynErr_KS2));max(abs(Res.dynErr_KS3))];
end
if case_nr == 1
    Res.max_dynErr_KS = [Res.max_dynErr_KS;max(abs(Res.dynErr_KS4))];
end
Res.max_dynErr_R = max(abs(Res.dynErr_R));
Res.max_dynErr_DR = max(abs(Res.dynErr_DR));
Res.max_dynErr_PPA = [max(abs(Res.dynErr_PPA0));max(abs(Res.dynErr_PPA1));max(abs(Res.dynErr_PPA2));max(abs(Res.dynErr_PPA3))];
Res.max_dynErr_PFI = [max(abs(Res.dynErr_PFI0));max(abs(Res.dynErr_PFI1));max(abs(Res.dynErr_PFI2));max(abs(Res.dynErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.max_dynErr_PPA = [Res.max_dynErr_PPA;max(abs(Res.dynErr_PPA4))];
    Res.max_dynErr_PFI = [Res.max_dynErr_PFI;max(abs(Res.dynErr_PFI4))];
end

Res.LoMErr_KS1 = (Res.K_pr_LoM_KS(1,1:end-1)-Res.K_KS(1,2:end))'./Res.K_KS(1,2:end)';
if case_nr<3
    Res.LoMErr_KS2 = (Res.K_pr_LoM_KS(2,1:end-1)-Res.K_KS(2,2:end))'./Res.K_KS(2,2:end)';
    Res.LoMErr_KS3 = (Res.K_pr_LoM_KS(3,1:end-1)-Res.K_KS(3,2:end))'./Res.K_KS(3,2:end)';
end
if case_nr == 1
    Res.LoMErr_KS4 = (Res.K_pr_LoM_KS(4,1:end-1)-Res.K_KS(4,2:end))'./Res.K_KS(4,2:end)';
end
Res.LoMErr_DR = (Res.K_pr_LoM_DR(1:end-1)-Res.K_DR(2:end))./Res.K_DR(2:end);
Res.LoMErr_PPA0 = squeeze((Res.K_pr_LoM(1,1,1:end-1)-Res.K(1,1,2:end))./Res.K(1,1,2:end));
Res.LoMErr_PPA1 = squeeze((Res.K_pr_LoM(1,2,1:end-1)-Res.K(1,2,2:end))./Res.K(1,2,2:end));
Res.LoMErr_PPA2 = squeeze((Res.K_pr_LoM(1,3,1:end-1)-Res.K(1,3,2:end))./Res.K(1,3,2:end));
Res.LoMErr_PPA3 = squeeze((Res.K_pr_LoM(1,4,1:end-1)-Res.K(1,4,2:end))./Res.K(1,4,2:end));
if case_nr==2 || case_nr>3
    Res.LoMErr_PPA4 = squeeze((Res.K_pr_LoM(1,5,1:end-1)-Res.K(1,5,2:end))./Res.K(1,5,2:end));
end
Res.LoMErr_PFI0 = squeeze((Res.K_pr_LoM(2,1,1:end-1)-Res.K(2,1,2:end))./Res.K(2,1,2:end));
Res.LoMErr_PFI1 = squeeze((Res.K_pr_LoM(2,2,1:end-1)-Res.K(2,2,2:end))./Res.K(2,2,2:end));
Res.LoMErr_PFI2 = squeeze((Res.K_pr_LoM(2,3,1:end-1)-Res.K(2,3,2:end))./Res.K(2,3,2:end));
Res.LoMErr_PFI3 = squeeze((Res.K_pr_LoM(2,4,1:end-1)-Res.K(2,4,2:end))./Res.K(2,4,2:end));
if case_nr==2 || case_nr>3
    Res.LoMErr_PFI4 = squeeze((Res.K_pr_LoM(2,5,1:end-1)-Res.K(2,5,2:end))./Res.K(2,5,2:end));
end

Res.median_LoMErr_KS = median((Res.LoMErr_KS1));
if case_nr<3
    Res.median_LoMErr_KS = [Res.median_LoMErr_KS;median(abs(Res.LoMErr_KS2));median(abs(Res.LoMErr_KS3))];
end
if case_nr == 1
    Res.median_LoMErr_KS = [Res.median_LoMErr_KS;median(abs(Res.LoMErr_KS4))];
end
Res.median_LoMErr_DR = median((Res.LoMErr_DR));
Res.median_LoMErr_PPA = [median((Res.LoMErr_PPA0));median((Res.LoMErr_PPA1));median((Res.LoMErr_PPA2));median((Res.LoMErr_PPA3))];
Res.median_LoMErr_PFI = [median((Res.LoMErr_PFI0));median((Res.LoMErr_PFI1));median((Res.LoMErr_PFI2));median((Res.LoMErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.median_LoMErr_PPA = [Res.median_LoMErr_PPA;median((Res.LoMErr_PPA4))];
    Res.median_LoMErr_PFI = [Res.median_LoMErr_PFI;median((Res.LoMErr_PFI4))];
end
Res.mean_LoMErr_KS = mean((Res.LoMErr_KS1));
if case_nr<3
    Res.mean_LoMErr_KS = [Res.mean_LoMErr_KS;mean(abs(Res.LoMErr_KS2));mean(abs(Res.LoMErr_KS3))];
end
if case_nr == 1
    Res.mean_LoMErr_KS = [Res.mean_LoMErr_KS;mean(abs(Res.LoMErr_KS4))];
end
Res.mean_LoMErr_DR = mean((Res.LoMErr_DR));
Res.mean_LoMErr_PPA = [mean((Res.LoMErr_PPA0));mean((Res.LoMErr_PPA1));mean((Res.LoMErr_PPA2));mean((Res.LoMErr_PPA3))];
Res.mean_LoMErr_PFI = [mean((Res.LoMErr_PFI0));mean((Res.LoMErr_PFI1));mean((Res.LoMErr_PFI2));mean((Res.LoMErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.mean_LoMErr_PPA = [Res.mean_LoMErr_PPA;mean((Res.LoMErr_PPA4))];
    Res.mean_LoMErr_PFI = [Res.mean_LoMErr_PFI;mean((Res.LoMErr_PFI4))];
end
Res.min_LoMErr_KS = min((Res.LoMErr_KS1));
if case_nr<3
    Res.min_LoMErr_KS = [Res.min_LoMErr_KS;min(abs(Res.LoMErr_KS2));min(abs(Res.LoMErr_KS3))];
end
if case_nr == 1
    Res.min_LoMErr_KS = [Res.min_LoMErr_KS;min(abs(Res.LoMErr_KS4))];
end
Res.min_LoMErr_DR = min((Res.LoMErr_DR));
Res.min_LoMErr_PPA = [min((Res.LoMErr_PPA0));min((Res.LoMErr_PPA1));min((Res.LoMErr_PPA2));min((Res.LoMErr_PPA3))];
Res.min_LoMErr_PFI = [min((Res.LoMErr_PFI0));min((Res.LoMErr_PFI1));min((Res.LoMErr_PFI2));min((Res.LoMErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.min_LoMErr_PPA = [Res.min_LoMErr_PPA;min((Res.LoMErr_PPA4))];
    Res.min_LoMErr_PFI = [Res.min_LoMErr_PFI;min((Res.LoMErr_PFI4))];
end
Res.max_LoMErr_KS = max((Res.LoMErr_KS1));
if case_nr<3
    Res.max_LoMErr_KS = [Res.max_LoMErr_KS;max(abs(Res.LoMErr_KS2));max(abs(Res.LoMErr_KS3))];
end
if case_nr == 1
    Res.max_LoMErr_KS = [Res.max_LoMErr_KS;max(abs(Res.LoMErr_KS4))];
end
Res.max_LoMErr_DR = max((Res.LoMErr_DR));
Res.max_LoMErr_PPA = [max((Res.LoMErr_PPA0));max((Res.LoMErr_PPA1));max((Res.LoMErr_PPA2));max((Res.LoMErr_PPA3))];
Res.max_LoMErr_PFI = [max((Res.LoMErr_PFI0));max((Res.LoMErr_PFI1));max((Res.LoMErr_PFI2));max((Res.LoMErr_PFI3))];
if case_nr==2 || case_nr>3
    Res.max_LoMErr_PPA = [Res.max_LoMErr_PPA;max((Res.LoMErr_PPA4))];
    Res.max_LoMErr_PFI = [Res.max_LoMErr_PFI;max((Res.LoMErr_PFI4))];
end

save(strcat(folder,'/Res.mat'),'-struct','Res');
end