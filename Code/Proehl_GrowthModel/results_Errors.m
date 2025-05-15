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
function results_Errors(case_nr,StaticParams,folder)

% truncation order 
M = 3;
M_KS = 4;

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------      
% Simulate the aggregate shock_____________________________________________
% Keep this section commented out to reproduce the results from the paper!

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
y=[-5e-4 0 0 -5e-4];
pdf_plot = 2*StaticParams.p(1:2)'*pdf_start;
pdf_plot(2:end) = pdf_plot(2:end)./(StaticParams.kGrid(2:end)-StaticParams.kGrid(1:end-1));
figure(case_nr)
p4=fill(x,y,'g','EdgeColor','none');
set(gca,'nextplot','add');
p5=fill(x,y,'r','EdgeColor','none');
set(gca,'nextplot','add');
p1=plot(StaticParams.kGrid,pdf_plot,'m','LineWidth',2);
set(gca,'nextplot','add');
p3=plot(StaticParams.kGrid,pdf_plot,'b','LineWidth',2);
set(gca,'nextplot','replacechildren');
legend([p1,p3,p4,p5],'Krusell-Smith', 'Proximal Point Algo', 'boom', 'recession','Location','best','AutoUpdate','off')
legend('boxoff')
axis([0 k_max -5e-4 3.5e-2])
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
K_KS = zeros(M_KS,Res.nsamples);
sd_KS = zeros(M_KS,Res.nsamples);
skew_KS = zeros(M_KS,Res.nsamples);
kurt_KS = zeros(M_KS,Res.nsamples);
gini_KS = zeros(M_KS,Res.nsamples);
dyn_euler_K_KS = zeros(M_KS,Res.nsamples);
K_pr_LoM_KS = zeros(M_KS,Res.nsamples);
idx_posProb_KS = zeros(M_KS,2,length(Res.kGrid));
idx = max(find(squeeze(pdf_KS(1,1,:))>eps,1,'last'),find(squeeze(pdf_KS(1,2,:))>eps,1,'last'));
idx_posProb_KS(:,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
dyn_idx_posProb_KS = idx_posProb_KS;

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
kb = zeros(2,M+1,Res.nsamples);
kg = zeros(2,M+1,Res.nsamples);
cb = zeros(2,M+1,Res.nsamples);
cg = zeros(2,M+1,Res.nsamples);
K = zeros(2,M+1,Res.nsamples);
sd = zeros(2,M+1,Res.nsamples);
skew = zeros(2,M+1,Res.nsamples);
kurt = zeros(2,M+1,Res.nsamples);
gini = zeros(2,M+1,Res.nsamples);
K_pr_LoM = zeros(2,M+1,Res.nsamples);
dyn_euler_K = zeros(2,M+1,Res.nsamples);
idx_posProb = zeros(2,M+1,2,length(Res.kGrid));
idx = max(find(squeeze(pdf(1,1,1,:))>eps,1,'last'),find(squeeze(pdf(1,1,2,:))>eps,1,'last'));
idx_posProb(:,:,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
dyn_idx_posProb = idx_posProb;

Res.stErr_KS1 = [];
Res.stErr_KS2 = [];
Res.stErr_KS3 = [];
Res.stErr_KS4 = [];
Res.stErr_PPA0 = [];
Res.stErr_PPA1 = [];
Res.stErr_PPA2 = [];
Res.stErr_PPA3 = [];
if M>3
    Res.stErr_PPA4 = [];
end
Res.stErr_PFI0 = [];
Res.stErr_PFI1 = [];
Res.stErr_PFI2 = [];
Res.stErr_PFI3 = [];
if M>3
    Res.stErr_PFI4 = [];
end

Res.dynErr_KS1 = [];
Res.dynErr_KS2 = [];
Res.dynErr_KS3 = [];
Res.dynErr_KS4 = [];
Res.dynErr_PPA0 = [];
Res.dynErr_PPA1 = [];
Res.dynErr_PPA2 = [];
Res.dynErr_PPA3 = [];
if M>3
    Res.dynErr_PPA4 = [];
end
Res.dynErr_PFI0 = [];
Res.dynErr_PFI1 = [];
Res.dynErr_PFI2 = [];
Res.dynErr_PFI3 = [];
if M>3
    Res.dynErr_PFI4 = [];
end

Res.stErr_KS1_algoGrid = [];
Res.stErr_KS2_algoGrid = [];
Res.stErr_KS3_algoGrid = [];
Res.stErr_KS4_algoGrid = [];
Res.stErr_PPA0_algoGrid = [];
Res.stErr_PPA1_algoGrid = [];
Res.stErr_PPA2_algoGrid = [];
Res.stErr_PPA3_algoGrid = [];
if M>3
    Res.stErr_PPA4_algoGrid = [];
end
Res.stErr_PFI0_algoGrid = [];
Res.stErr_PFI1_algoGrid = [];
Res.stErr_PFI2_algoGrid = [];
Res.stErr_PFI3_algoGrid = [];
if M>3
    Res.stErr_PFI4_algoGrid = [];
end

Res.dynErr_KS1_algoGrid = [];
Res.dynErr_KS2_algoGrid = [];
Res.dynErr_KS3_algoGrid = [];
Res.dynErr_KS4_algoGrid = [];
Res.dynErr_PPA0_algoGrid = [];
Res.dynErr_PPA1_algoGrid = [];
Res.dynErr_PPA2_algoGrid = [];
Res.dynErr_PPA3_algoGrid = [];
if M>3
    Res.dynErr_PPA4_algoGrid = [];
end
Res.dynErr_PFI0_algoGrid = [];
Res.dynErr_PFI1_algoGrid = [];
Res.dynErr_PFI2_algoGrid = [];
Res.dynErr_PFI3_algoGrid = [];
if M>3
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
                  [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))')])',Res.kGrid)';
        c_prime = interp1(StaticParams.kGrid_pol,squeeze(...
                  [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))')])',Res.kGrid)';   
        if i==1
            kb_0 = K(method+1,order_ct+1,i);
            kg_0 = K(method+1,order_ct+1,i);
        else
            kb_0 = max(StaticParams.k_min,kb(method+1,order_ct+1,i-1));
            kg_0 = max(StaticParams.k_min,kg(method+1,order_ct+1,i-1));        
        end
        kb(method+1,order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))'))',kb_0)';
        kg(method+1,order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))'))',kg_0)';
        cb(method+1,order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))'))',kb_0)';   
        cg(method+1,order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct+1,1:(Poly.phi_N+1)))'))',kg_0)';   
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
        euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                      1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                      .*(interp1(StaticParams.kGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')')...
                      .^(-StaticParams.gamma);
        euler_c_prime = min(StaticParams.wealth(Res.z_ag(i),squeeze(K(method+1,order_ct+1,i)),Res.kGrid),...
                                    [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime(1,:))).^(-1/StaticParams.gamma);...
                                     (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime(2,:))).^(-1/StaticParams.gamma)]);
        % standard error
        Err = log10(abs(1-euler_c_prime./c_prime));%(c_prime-euler_c_prime)./max(eps,euler_c_prime);
        Err_algoGrid = log10(abs(1-euler_c_prime(:,algo_idx)./c_prime(:,algo_idx)));%(c_prime(:,algo_idx)-euler_c_prime(:,algo_idx))./max(eps,euler_c_prime(:,algo_idx));
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
        k_prime_next = [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(dyn_weights(method+1,order_ct+1,1:(Poly.phi_N+1)))');...
                        getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(dyn_weights(method+1,order_ct+1,1:(Poly.phi_N+1)))')];  
        euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                  *sum(interp1(StaticParams.kGrid_pol,sort(k_prime_next,2)',StaticParams.kGrid)'.*squeeze(dyn_pdf(method+1,order_ct+1,:,:)),2);
        euler_weights = calcWeights(squeeze(dyn_cdf(method+1,order_ct+1,:,:)),Poly,StaticParams,Res.z_ag(i),...
                        interp1(StaticParams.kGrid_pol,sort(k_prime_next,2)',StaticParams.kGrid)');
        euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                      1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                      .*(interp1(StaticParams.kGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')')...
                      .^(-StaticParams.gamma);
        dyn_euler_c_prime = min(StaticParams.wealth(Res.z_ag(i),squeeze(dyn_euler_K(method+1,order_ct+1,i)),Res.kGrid),...
                                    [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(interp1(StaticParams.kGrid_pol,sort(k_prime_next(1,:),2)',Res.kGrid)')).^(-1/StaticParams.gamma);...
                                     (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(interp1(StaticParams.kGrid_pol,sort(k_prime_next(2,:),2)',Res.kGrid)')).^(-1/StaticParams.gamma)]);
        dyn_euler_k_prime = StaticParams.wealth(Res.z_ag(i),squeeze(dyn_euler_K(method+1,order_ct+1,i)),Res.kGrid)-dyn_euler_c_prime;
        % dynamic error
        Err = log10(abs(1-dyn_euler_c_prime./c_prime));%(c_prime-dyn_euler_c_prime)./max(eps,dyn_euler_c_prime);
        Err_algoGrid = log10(abs(1-dyn_euler_c_prime(:,algo_idx)./c_prime(:,algo_idx)));%(c_prime(:,algo_idx)-dyn_euler_c_prime(:,algo_idx))./max(eps,dyn_euler_c_prime(:,algo_idx));
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
    for order_ct = 1:M_KS
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
        switch order_ct
            case 1
            k_prime_KS = max(0,[interpn(KS.k,KS.km,KS.kprime(:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))),'spline')';...
                             interpn(KS.k,KS.km,KS.kprime(:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))),'spline')']);
            case 2
            k_prime_KS = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_KS(order_ct,i))),'spline')';...
                             interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_KS(order_ct,i))),'spline')']);
            case 3
            k_prime_KS = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_KS(order_ct,i))),min(KS.kskew_max,max(KS.kskew_min,skew_KS(order_ct,i))),'spline')';...
                             interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_KS(order_ct,i))),min(KS.kskew_max,max(KS.kskew_min,skew_KS(order_ct,i))),'spline')']);
            case 4
            k_prime_KS = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_KS(order_ct,i))),min(KS.kskew_max,max(KS.kskew_min,skew_KS(order_ct,i))),min(KS.kkurt_max,max(KS.kkurt_min,kurt_KS(order_ct,i))),'spline')';...
                             interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_KS(order_ct,i))),min(KS.kskew_max,max(KS.kskew_min,skew_KS(order_ct,i))),min(KS.kkurt_max,max(KS.kkurt_min,kurt_KS(order_ct,i))),'spline')']);
        end
        c_prime_KS = max(0,StaticParams.wealth(Res.z_ag(i),K_KS(order_ct,i),Res.kGrid)-k_prime_KS);
        % standard Euler policies 
        switch order_ct
            case 1
            K_pr_LoM_KS(order_ct,i) = exp(KS.B1((order_ct+1)*Res.z_ag(i)+1)+KS.B1((order_ct+1)*Res.z_ag(i)+2)*log(K_KS(order_ct,i)));
            case 2
            K_pr_LoM_KS(order_ct,i) = exp(KS.B1((order_ct+1)*Res.z_ag(i)+1)+KS.B1((order_ct+1)*Res.z_ag(i)+2)*log(K_KS(order_ct,i))+KS.B1((order_ct+1)*Res.z_ag(i)+3)*log(sd_KS(order_ct,i)));
            case 3
            K_pr_LoM_KS(order_ct,i) = exp(KS.B1((order_ct+1)*Res.z_ag(i)+1)+KS.B1((order_ct+1)*Res.z_ag(i)+2)*log(K_KS(order_ct,i))+KS.B1((order_ct+1)*Res.z_ag(i)+3)*log(sd_KS(order_ct,i))+KS.B1((order_ct+1)*Res.z_ag(i)+4)*(skew_KS(order_ct,i)));
            case 4
            K_pr_LoM_KS(order_ct,i) = exp(KS.B1((order_ct+1)*Res.z_ag(i)+1)+KS.B1((order_ct+1)*Res.z_ag(i)+2)*log(K_KS(order_ct,i))+KS.B1((order_ct+1)*Res.z_ag(i)+3)*log(sd_KS(order_ct,i))+KS.B1((order_ct+1)*Res.z_ag(i)+4)*(skew_KS(order_ct,i))+KS.B1((order_ct+1)*Res.z_ag(i)+5)*log(kurt_KS(order_ct,i)));
        end
        euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                  *sum(interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'.*squeeze(pdf_KS(order_ct,:,:)),2);
        switch order_ct
            case 1
            euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                          1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                          .*max(1e-15,([StaticParams.wealth(0,euler_K,grid);StaticParams.wealth(1,euler_K,grid)]...
                                -interp1(StaticParams.kGrid,...
                                  [interpn(KS.k,KS.km,KS.kprime(:,:,1,1),StaticParams.kGrid,euler_K,'spline')';...
                                   interpn(KS.k,KS.km,KS.kprime(:,:,1,2),StaticParams.kGrid,euler_K,'spline')';...
                                   interpn(KS.k,KS.km,KS.kprime(:,:,2,1),StaticParams.kGrid,euler_K,'spline')';...
                                   interpn(KS.k,KS.km,KS.kprime(:,:,2,2),StaticParams.kGrid,euler_K,'spline')']',...
                                 grid,'linear','extrap')'))...
                          .^(-StaticParams.gamma);
            case 2        
            euler_sd = sqrt([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^2 ... 
                            .*squeeze(pdf_KS(order_ct,:,:)),2));
            euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                          1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                          .*max(1e-15,([StaticParams.wealth(0,euler_K,grid);StaticParams.wealth(1,euler_K,grid)]...
                                -interp1(StaticParams.kGrid,...
                                  [interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,1,1),StaticParams.kGrid,euler_K,euler_sd,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,1,2),StaticParams.kGrid,euler_K,euler_sd,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,2,1),StaticParams.kGrid,euler_K,euler_sd,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,2,2),StaticParams.kGrid,euler_K,euler_sd,'spline')']',...
                                 grid,'linear','extrap')'))...
                          .^(-StaticParams.gamma);
            case 3
            euler_sd = sqrt([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^2 ... 
                            .*squeeze(pdf_KS(order_ct,:,:)),2));
            euler_skew = ([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^3 ... 
                            .*squeeze(pdf_KS(order_ct,:,:)),2))./(euler_sd).^3;
            euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                          1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                          .*max(1e-15,([StaticParams.wealth(0,euler_K,grid);StaticParams.wealth(1,euler_K,grid)]...
                                -interp1(StaticParams.kGrid,...
                                  [interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,1,1),StaticParams.kGrid,euler_K,euler_sd,euler_skew,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,1,2),StaticParams.kGrid,euler_K,euler_sd,euler_skew,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,2,1),StaticParams.kGrid,euler_K,euler_sd,euler_skew,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,2,2),StaticParams.kGrid,euler_K,euler_sd,euler_skew,'spline')']',...
                                 grid,'linear','extrap')'))...
                          .^(-StaticParams.gamma);
            case 4
            euler_sd = sqrt([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^2 ... 
                            .*squeeze(pdf_KS(order_ct,:,:)),2));
            euler_skew = ([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^3 ... 
                            .*squeeze(pdf_KS(order_ct,:,:)),2))./(euler_sd).^3;
            euler_kurt = ([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^4 ... 
                            .*squeeze(pdf_KS(order_ct,:,:)),2))./(euler_sd).^4;
            euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                          1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                          .*max(1e-15,([StaticParams.wealth(0,euler_K,grid);StaticParams.wealth(1,euler_K,grid)]...
                                -interp1(StaticParams.kGrid,...
                                  [interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,1,1),StaticParams.kGrid,euler_K,euler_sd,euler_skew,euler_kurt,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,1,2),StaticParams.kGrid,euler_K,euler_sd,euler_skew,euler_kurt,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,2,1),StaticParams.kGrid,euler_K,euler_sd,euler_skew,euler_kurt,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,2,2),StaticParams.kGrid,euler_K,euler_sd,euler_skew,euler_kurt,'spline')']',...
                                 grid,'linear','extrap')'))...
                          .^(-StaticParams.gamma);
        end
        euler_c_prime_KS = min(StaticParams.wealth(Res.z_ag(i),K_KS(order_ct,i),Res.kGrid),...
                           [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime_KS(1,:))).^(-1/StaticParams.gamma);...
                            (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime_KS(2,:))).^(-1/StaticParams.gamma)]);
        % standard error
        Err = log10(abs(1-euler_c_prime_KS./c_prime_KS));%(c_prime_KS-euler_c_prime_KS)./max(eps,euler_c_prime_KS);
        Err_algoGrid = log10(abs(1-euler_c_prime_KS(:,algo_idx)./c_prime_KS(:,algo_idx)));%(c_prime_KS(:,algo_idx)-euler_c_prime_KS(:,algo_idx))./max(eps,euler_c_prime_KS(:,algo_idx));
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
            idx_posProb_KS(order_ct,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
            pdf_KS(order_ct,:,:) = pdf_KS(order_ct,:,:).*(pdf_KS(order_ct,:,:)>eps);
            pdf_KS(order_ct,:,:) = pdf_KS(order_ct,:,:)./repmat(sum(pdf_KS(order_ct,:,:),3),[1,1,size(pdf_KS,3)]);   
            cdf_KS(order_ct,:,:) = cummax(min(1,cumsum(pdf_KS(order_ct,:,:),3)),3);
            cdf_KS(order_ct,1,cdf_KS(order_ct,1,:)==max(cdf_KS(order_ct,1,:))) = 1;
            cdf_KS(order_ct,2,cdf_KS(order_ct,2,:)==max(cdf_KS(order_ct,2,:))) = 1; 
        end
        
        % dynamic Euler policies and distributions
        dyn_euler_K_KS(order_ct,i) = StaticParams.kGrid*squeeze(dyn_pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
        switch order_ct
            case 1
            k_prime_next = max(0,[interpn(KS.k,KS.km,KS.kprime(:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))),'spline')';...
                             interpn(KS.k,KS.km,KS.kprime(:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))),'spline')']);
            case 2
            sd_aux = sqrt(StaticParams.kGrid.^2*squeeze(dyn_pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]);
            k_prime_next = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_aux)),'spline')';...
                             interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_aux)),'spline')']);
            case 3
            sd_aux = sqrt(StaticParams.kGrid.^2*squeeze(dyn_pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]);
            skew_aux = ((StaticParams.kGrid-squeeze(K_KS(order_ct,i))).^3 ...
                         *squeeze(dyn_pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                         ./squeeze(sd_aux)^3;
            k_prime_next = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_aux)),min(KS.kskew_max,max(KS.kskew_min,skew_aux)),'spline')';...
                             interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_aux)),min(KS.kskew_max,max(KS.kskew_min,skew_aux)),'spline')']);
            case 4
            sd_aux = sqrt(StaticParams.kGrid.^2*squeeze(dyn_pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]);
            skew_aux = ((StaticParams.kGrid-squeeze(K_KS(order_ct,i))).^3 ...
                         *squeeze(dyn_pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                         ./squeeze(sd_aux)^3;
            kurt_aux = ((StaticParams.kGrid-squeeze(K_KS(order_ct,i))).^4 ...
                         *squeeze(dyn_pdf_KS(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                         ./squeeze(sd_aux)^4;
            k_prime_next = max(0,[interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,Res.z_ag(i)+1,1),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_aux)),min(KS.kskew_max,max(KS.kskew_min,skew_aux)),min(KS.kkurt_max,max(KS.kkurt_min,kurt_aux)),'spline')';...
                             interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,Res.z_ag(i)+1,2),Res.kGrid,min(KS.km_max,max(KS.km_min,dyn_euler_K_KS(order_ct,i))),min(KS.kvar_max,max(KS.kvar_min,sd_aux)),min(KS.kskew_max,max(KS.kskew_min,skew_aux)),min(KS.kkurt_max,max(KS.kkurt_min,kurt_aux)),'spline')']);
        end
        euler_K = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                  *sum(interp1(Res.kGrid,sort(k_prime_next,2)',StaticParams.kGrid)'.*squeeze(dyn_pdf_KS(order_ct,:,:)),2);
        switch order_ct
            case 1
            euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                          1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                          .*max(1e-15,([StaticParams.wealth(0,euler_K,grid);StaticParams.wealth(1,euler_K,grid)]...
                                -interp1(StaticParams.kGrid,...
                                  [interpn(KS.k,KS.km,KS.kprime(:,:,1,1),StaticParams.kGrid,euler_K,'spline')';...
                                   interpn(KS.k,KS.km,KS.kprime(:,:,1,2),StaticParams.kGrid,euler_K,'spline')';...
                                   interpn(KS.k,KS.km,KS.kprime(:,:,2,1),StaticParams.kGrid,euler_K,'spline')';...
                                   interpn(KS.k,KS.km,KS.kprime(:,:,2,2),StaticParams.kGrid,euler_K,'spline')']',...
                                 grid,'linear','extrap')'))...
                          .^(-StaticParams.gamma);
            case 2        
            euler_sd = sqrt([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^2 ... 
                            .*squeeze(dyn_pdf_KS(order_ct,:,:)),2));
            euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                          1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                          .*max(1e-15,([StaticParams.wealth(0,euler_K,grid);StaticParams.wealth(1,euler_K,grid)]...
                                -interp1(StaticParams.kGrid,...
                                  [interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,1,1),StaticParams.kGrid,euler_K,euler_sd,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,1,2),StaticParams.kGrid,euler_K,euler_sd,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,2,1),StaticParams.kGrid,euler_K,euler_sd,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kprime(:,:,:,2,2),StaticParams.kGrid,euler_K,euler_sd,'spline')']',...
                                 grid,'linear','extrap')'))...
                          .^(-StaticParams.gamma);
            case 3
            euler_sd = sqrt([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^2 ... 
                            .*squeeze(dyn_pdf_KS(order_ct,:,:)),2));
            euler_skew = ([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^3 ... 
                            .*squeeze(dyn_pdf_KS(order_ct,:,:)),2))./(euler_sd).^3;
            euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                          1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                          .*max(1e-15,([StaticParams.wealth(0,euler_K,grid);StaticParams.wealth(1,euler_K,grid)]...
                                -interp1(StaticParams.kGrid,...
                                  [interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,1,1),StaticParams.kGrid,euler_K,euler_sd,euler_skew,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,1,2),StaticParams.kGrid,euler_K,euler_sd,euler_skew,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,2,1),StaticParams.kGrid,euler_K,euler_sd,euler_skew,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kprime(:,:,:,:,2,2),StaticParams.kGrid,euler_K,euler_sd,euler_skew,'spline')']',...
                                 grid,'linear','extrap')'))...
                          .^(-StaticParams.gamma);
            case 4
            euler_sd = sqrt([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^2 ... 
                            .*squeeze(dyn_pdf_KS(order_ct,:,:)),2));
            euler_skew = ([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^3 ... 
                            .*squeeze(dyn_pdf_KS(order_ct,:,:)),2))./(euler_sd).^3;
            euler_kurt = ([StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]...
                       *sum((interp1(Res.kGrid,sort(k_prime_KS,2)',StaticParams.kGrid)'-euler_K).^4 ... 
                            .*squeeze(dyn_pdf_KS(order_ct,:,:)),2))./(euler_sd).^4;
            euler_expec = @(grid) repmat([1-StaticParams.delta+StaticParams.rate(0,euler_K);...
                                          1-StaticParams.delta+StaticParams.rate(1,euler_K)],[1,length(grid)])...
                          .*max(1e-15,([StaticParams.wealth(0,euler_K,grid);StaticParams.wealth(1,euler_K,grid)]...
                                -interp1(StaticParams.kGrid,...
                                  [interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,1,1),StaticParams.kGrid,euler_K,euler_sd,euler_skew,euler_kurt,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,1,2),StaticParams.kGrid,euler_K,euler_sd,euler_skew,euler_kurt,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,2,1),StaticParams.kGrid,euler_K,euler_sd,euler_skew,euler_kurt,'spline')';...
                                   interpn(KS.k,KS.km,KS.kvar,KS.kskew,KS.kkurt,KS.kprime(:,:,:,:,:,2,2),StaticParams.kGrid,euler_K,euler_sd,euler_skew,euler_kurt,'spline')']',...
                                 grid,'linear','extrap')'))...
                          .^(-StaticParams.gamma);
        end
        dyn_euler_c_prime_KS = min(StaticParams.wealth(Res.z_ag(i),dyn_euler_K_KS(i),Res.kGrid),...
                               [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime_next(1,:))).^(-1/StaticParams.gamma);...
                                (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime_next(2,:))).^(-1/StaticParams.gamma)]);
        dyn_euler_k_prime_KS = StaticParams.wealth(Res.z_ag(i),dyn_euler_K_KS(order_ct,i),Res.kGrid)-dyn_euler_c_prime_KS;
        % dynamic error
        Err = log10(abs(1-dyn_euler_c_prime_KS./c_prime_KS));%(c_prime_KS-dyn_euler_c_prime_KS)./max(eps,dyn_euler_c_prime_KS);
        Err_algoGrid = log10(abs(1-dyn_euler_c_prime_KS(:,algo_idx)./c_prime_KS(:,algo_idx)));%(c_prime_KS(:,algo_idx)-dyn_euler_c_prime_KS(:,algo_idx))./max(eps,dyn_euler_c_prime_KS(:,algo_idx));
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
            dyn_idx_posProb_KS(order_ct,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
            dyn_pdf_KS(order_ct,:,:) = dyn_pdf_KS(order_ct,:,:).*(dyn_pdf_KS(order_ct,:,:)>eps);
            dyn_pdf_KS(order_ct,:,:) = dyn_pdf_KS(order_ct,:,:)./repmat(sum(dyn_pdf_KS(order_ct,:,:),3),[1,1,size(dyn_pdf_KS,3)]);
            dyn_cdf_KS(order_ct,:,:) = cummax(min(1,cumsum(dyn_pdf_KS(order_ct,:,:),3)),3);
            dyn_cdf_KS(order_ct,1,dyn_cdf_KS(order_ct,1,:)==max(dyn_cdf_KS(order_ct,1,:))) = 1;
            dyn_cdf_KS(order_ct,2,dyn_cdf_KS(order_ct,2,:)==max(dyn_cdf_KS(order_ct,2,:))) = 1;  
        end
    end
    % update video
    figure(case_nr)
    if Res.z_ag(i)==1
        p5=fill(x,y,'r','EdgeColor','none');
        set(gca,'nextplot','add');
        p4=fill(x,y,'g','EdgeColor','none');
    else
        p4=fill(x,y,'g','EdgeColor','none');
        set(gca,'nextplot','add');
        p5=fill(x,y,'r','EdgeColor','none');
    end
    set(gca,'nextplot','add');
    pdf_plot = 2*StaticParams.p(Res.z_ag(i)*2+(1:2))'*squeeze(pdf_KS(1,:,:));
    pdf_plot(2:end) = pdf_plot(2:end)./(StaticParams.kGrid(2:end)-StaticParams.kGrid(1:end-1));
    p1=plot(StaticParams.kGrid,pdf_plot,'m','LineWidth',2);
    set(gca,'nextplot','add');
    pdf_plot = 2*StaticParams.p(Res.z_ag(i)*2+(1:2))'*squeeze(pdf(1,end,:,:));
    pdf_plot(2:end) = pdf_plot(2:end)./(StaticParams.kGrid(2:end)-StaticParams.kGrid(1:end-1));
    p3=plot(StaticParams.kGrid,pdf_plot,'b','LineWidth',2);
    set(gca,'nextplot','replacechildren');
    legend([p1,p3,p4,p5],'Krusell-Smith', 'Proximal Point Algo', 'boom', 'recession','Location','best','AutoUpdate','off')
    legend('boxoff')
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

Res.kb = kb;
Res.kg = kg;
Res.cb = cb;
Res.cg = cg;

Res.K_KS = K_KS;
Res.sd_KS = sd_KS;
Res.skew_KS = skew_KS;
Res.kurt_KS = kurt_KS;
Res.gini_KS = gini_KS;
Res.K_pr_LoM_KS = K_pr_LoM_KS;
Res.dyn_euler_K_KS = dyn_euler_K_KS;

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
Res.median_stErr_KS = [];
Res.median_stErr_PPA = [];
Res.median_stErr_PFI = [];
Res.mean_stErr_KS = [];
Res.mean_stErr_PPA = [];
Res.mean_stErr_PFI = [];
Res.min_stErr_KS = [];
Res.min_stErr_PPA = [];
Res.min_stErr_PFI = [];
Res.max_stErr_KS = [];
Res.max_stErr_PPA = [];
Res.max_stErr_PFI = [];

Res.median_dynErr_KS = [];
Res.median_dynErr_PPA = [];
Res.median_dynErr_PFI = [];
Res.mean_dynErr_KS = [];
Res.mean_dynErr_PPA = [];
Res.mean_dynErr_PFI = [];
Res.min_dynErr_KS = [];
Res.min_dynErr_PPA = [];
Res.min_dynErr_PFI = [];
Res.max_dynErr_KS = [];
Res.max_dynErr_PPA = [];
Res.max_dynErr_PFI = [];
for i=0:M
    Res.median_stErr_KS = [Res.median_stErr_KS;median(eval(strcat('Res.stErr_KS',num2str(i+1))))];
    Res.median_stErr_PPA = [Res.median_stErr_PPA;median(eval(strcat('Res.stErr_PPA',num2str(i))))];
    Res.median_stErr_PFI = [Res.median_stErr_PFI;median(eval(strcat('Res.stErr_PFI',num2str(i))))];
    Res.mean_stErr_KS = [Res.mean_stErr_KS;mean(eval(strcat('Res.stErr_KS',num2str(i+1))))];
    Res.mean_stErr_PPA = [Res.mean_stErr_PPA;mean(eval(strcat('Res.stErr_PPA',num2str(i))))];
    Res.mean_stErr_PFI = [Res.mean_stErr_PFI;mean(eval(strcat('Res.stErr_PFI',num2str(i))))];
    Res.min_stErr_KS = [Res.min_stErr_KS;min(eval(strcat('Res.stErr_KS',num2str(i+1))))];
    Res.min_stErr_PPA = [Res.min_stErr_PPA;min(eval(strcat('Res.stErr_PPA',num2str(i))))];
    Res.min_stErr_PFI = [Res.min_stErr_PFI;min(eval(strcat('Res.stErr_PFI',num2str(i))))];
    Res.max_stErr_KS = [Res.max_stErr_KS;max(eval(strcat('Res.stErr_KS',num2str(i+1))))];
    Res.max_stErr_PPA = [Res.max_stErr_PPA;max(eval(strcat('Res.stErr_PPA',num2str(i))))];
    Res.max_stErr_PFI = [Res.max_stErr_PFI;max(eval(strcat('Res.stErr_PFI',num2str(i))))];
    
    Res.median_dynErr_KS = [Res.median_dynErr_KS;median(eval(strcat('Res.dynErr_KS',num2str(i+1))))];
    Res.median_dynErr_PPA = [Res.median_dynErr_PPA;median(eval(strcat('Res.dynErr_PPA',num2str(i))))];
    Res.median_dynErr_PFI = [Res.median_dynErr_PFI;median(eval(strcat('Res.dynErr_PFI',num2str(i))))];
    Res.mean_dynErr_KS = [Res.mean_dynErr_KS;mean(eval(strcat('Res.dynErr_KS',num2str(i+1))))];
    Res.mean_dynErr_PPA = [Res.mean_dynErr_PPA;mean(eval(strcat('Res.dynErr_PPA',num2str(i))))];
    Res.mean_dynErr_PFI = [Res.mean_dynErr_PFI;mean(eval(strcat('Res.dynErr_PFI',num2str(i))))];
    Res.min_dynErr_KS = [Res.min_dynErr_KS;min(eval(strcat('Res.dynErr_KS',num2str(i+1))))];
    Res.min_dynErr_PPA = [Res.min_dynErr_PPA;min(eval(strcat('Res.dynErr_PPA',num2str(i))))];
    Res.min_dynErr_PFI = [Res.min_dynErr_PFI;min(eval(strcat('Res.dynErr_PFI',num2str(i))))];
    Res.max_dynErr_KS = [Res.max_dynErr_KS;max(eval(strcat('Res.dynErr_KS',num2str(i+1))))];
    Res.max_dynErr_PPA = [Res.max_dynErr_PPA;max(eval(strcat('Res.dynErr_PPA',num2str(i))))];
    Res.max_dynErr_PFI = [Res.max_dynErr_PFI;max(eval(strcat('Res.dynErr_PFI',num2str(i))))];
end 

Res.LoMErr_KS1 = squeeze(log10(abs(Res.K_pr_LoM_KS(1,1:end-1)./Res.K_KS(1,2:end)-1)))';
Res.LoMErr_KS2 = squeeze(log10(abs(Res.K_pr_LoM_KS(2,1:end-1)./Res.K_KS(2,2:end)-1)))';
Res.LoMErr_KS3 = squeeze(log10(abs(Res.K_pr_LoM_KS(3,1:end-1)./Res.K_KS(3,2:end)-1)))';
Res.LoMErr_KS4 = squeeze(log10(abs(Res.K_pr_LoM_KS(4,1:end-1)./Res.K_KS(4,2:end)-1)))';
Res.LoMErr_PPA0 = squeeze(log10(abs(Res.K_pr_LoM(1,1,1:end-1)./Res.K(1,1,2:end)-1)));
Res.LoMErr_PPA1 = squeeze(log10(abs(Res.K_pr_LoM(1,2,1:end-1)./Res.K(1,2,2:end)-1)));
Res.LoMErr_PPA2 = squeeze(log10(abs(Res.K_pr_LoM(1,3,1:end-1)./Res.K(1,3,2:end)-1)));
Res.LoMErr_PPA3 = squeeze(log10(abs(Res.K_pr_LoM(1,4,1:end-1)./Res.K(1,4,2:end)-1)));
if M>3
    Res.LoMErr_PPA4 = squeeze(log10(abs(Res.K_pr_LoM(1,5,1:end-1)./Res.K(1,5,2:end)-1)));
end
Res.LoMErr_PFI0 = squeeze(log10(abs(Res.K_pr_LoM(2,1,1:end-1)./Res.K(2,1,2:end)-1)));
Res.LoMErr_PFI1 = squeeze(log10(abs(Res.K_pr_LoM(2,2,1:end-1)./Res.K(2,2,2:end)-1)));
Res.LoMErr_PFI2 = squeeze(log10(abs(Res.K_pr_LoM(2,3,1:end-1)./Res.K(2,3,2:end)-1)));
Res.LoMErr_PFI3 = squeeze(log10(abs(Res.K_pr_LoM(2,4,1:end-1)./Res.K(2,4,2:end)-1)));
if M>3
    Res.LoMErr_PFI4 = squeeze(log10(abs(Res.K_pr_LoM(2,5,1:end-1)./Res.K(2,5,2:end)-1)));
end

Res.median_LoMErr_KS = [];
Res.median_LoMErr_PPA = [];
Res.median_LoMErr_PFI = [];
Res.mean_LoMErr_KS = [];
Res.mean_LoMErr_PPA = [];
Res.mean_LoMErr_PFI = [];
Res.min_LoMErr_KS = [];
Res.min_LoMErr_PPA = [];
Res.min_LoMErr_PFI = [];
Res.max_LoMErr_KS = [];
Res.max_LoMErr_PPA = [];
Res.max_LoMErr_PFI = [];
for i=0:M
    Res.median_LoMErr_KS = [Res.median_LoMErr_KS;median(eval(strcat('Res.LoMErr_KS',num2str(i+1))))];
    Res.median_LoMErr_PPA = [Res.median_LoMErr_PPA;median(eval(strcat('Res.LoMErr_PPA',num2str(i))))];
    Res.median_LoMErr_PFI = [Res.median_LoMErr_PFI;median(eval(strcat('Res.LoMErr_PFI',num2str(i))))];
    Res.mean_LoMErr_KS = [Res.mean_LoMErr_KS;mean(eval(strcat('Res.LoMErr_KS',num2str(i+1))))];
    Res.mean_LoMErr_PPA = [Res.mean_LoMErr_PPA;mean(eval(strcat('Res.LoMErr_PPA',num2str(i))))];
    Res.mean_LoMErr_PFI = [Res.mean_LoMErr_PFI;mean(eval(strcat('Res.LoMErr_PFI',num2str(i))))];
    Res.min_LoMErr_KS = [Res.min_LoMErr_KS;min(eval(strcat('Res.LoMErr_KS',num2str(i+1))))];
    Res.min_LoMErr_PPA = [Res.min_LoMErr_PPA;min(eval(strcat('Res.LoMErr_PPA',num2str(i))))];
    Res.min_LoMErr_PFI = [Res.min_LoMErr_PFI;min(eval(strcat('Res.LoMErr_PFI',num2str(i))))];
    Res.max_LoMErr_KS = [Res.max_LoMErr_KS;max(eval(strcat('Res.LoMErr_KS',num2str(i+1))))];
    Res.max_LoMErr_PPA = [Res.max_LoMErr_PPA;max(eval(strcat('Res.LoMErr_PPA',num2str(i))))];
    Res.max_LoMErr_PFI = [Res.max_LoMErr_PFI;max(eval(strcat('Res.LoMErr_PFI',num2str(i))))];
end

save(strcat(folder,'/Res.mat'),'-struct','Res');
end