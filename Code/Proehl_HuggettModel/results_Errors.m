% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This file computes the expected ergodic distributions of the different 
% algorithm solutions.
%__________________________________________________________________________
function results_Errors(StaticParams,folder)

% truncation order 
M = 4;

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
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
Res.kGrid = unique([StaticParams.kGrid_pol,linspace(StaticParams.k_min,StaticParams.k_max,1001)]);
algo_idx = logical(sum(repmat(Res.kGrid,[length(StaticParams.kGrid_pol),1])==repmat(StaticParams.kGrid_pol',[1,length(Res.kGrid)]),1));
preSol = load(strcat(folder,'/preSol_PPA.mat'));
pdf_start = [[0.5,0.5]*[preSol.pdf_cond_b(1,:);preSol.pdf_cond_g(1,:)];...
             [0.5,0.5]*[preSol.pdf_cond_b(2,:);preSol.pdf_cond_g(2,:)]];
cdf_start = cummax(min(1,cumsum(pdf_start,2)),2);
cdf_start(1,cdf_start(1,:)==max(cdf_start(1,:)))=1;
cdf_start(2,cdf_start(2,:)==max(cdf_start(2,:)))=1;

v = VideoWriter(strcat(folder,'/Distr_Sim.avi'));
open(v);
k_max=0.75*StaticParams.k_max;
x=[StaticParams.k_min StaticParams.k_min k_max k_max];
y=[-0.015 0 0 -0.015];
[grid,idx_pdf] = sort([StaticParams.kGrid,StaticParams.kGrid(1:end-1)+1e-5]);
figure
p4=fill(x,y,'g','EdgeColor','none');
set(gca,'nextplot','add');
p5=fill(x,y,'r','EdgeColor','none');
set(gca,'nextplot','add');
pdf_aux = [pdf_start(1,:),pdf_start(1,2:end)];
pdf_aux = pdf_aux(idx_pdf);
p1=plot(grid,pdf_aux,'m','LineWidth',2);
set(gca,'nextplot','add');
pdf_aux = [pdf_start(2,:),pdf_start(2,2:end)];
pdf_aux = pdf_aux(idx_pdf);
p2=plot(grid,pdf_aux,'b','LineWidth',2);
set(gca,'nextplot','replacechildren');
legend([p1,p2,p4,p5],'low income', 'high income', 'boom', 'recession','Location','best','AutoUpdate','off')
legend('boxoff')
axis([StaticParams.k_min k_max -0.015 0.8])
title('Simulation of Asset Holding Distributions')
xlabel('a')
ylabel('probability')
frame = getframe(gcf);
writeVideo(v,frame);

% Set variables for PPA order 1 to M
pdf = repmat(permute(pdf_start,[4,3,1,2]),[2,M,1,1]);
cdf = repmat(permute(cdf_start,[4,3,1,2]),[2,M,1,1]);
weights = zeros(2,M,M);
for method = 0:1
for order_ct=1:M
    if method==0
        Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PPA.mat'));
    else
        Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
    end
    weights(method+1,order_ct,1:(Poly.phi_N+1))...
    = calcWeights(squeeze(cdf(method+1,order_ct,:,:)),Poly,StaticParams,0);
end
end
dyn_pdf = pdf;
dyn_cdf = cdf;
dyn_weights = weights;
kb = zeros(2,M,Res.nsamples,1);
kg = zeros(2,M,Res.nsamples,1);
cb = zeros(2,M,Res.nsamples,1);
cg = zeros(2,M,Res.nsamples,1);
pr = zeros(2,M,Res.nsamples,1);
K = zeros(2,M,Res.nsamples,1);
sd = zeros(2,M,Res.nsamples,1);
skew = zeros(2,M,Res.nsamples,1);
kurt = zeros(2,M,Res.nsamples,1);
dyn_euler_K = zeros(2,M,Res.nsamples,1);
idx_posProb = zeros(2,M,2,length(Res.kGrid));
idx = max(find(pdf(1,1,1,:)>eps,1,'last'),find(pdf(1,1,2,:)>eps,1,'last'));
idx_posProb(:,:,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
dyn_idx_posProb = idx_posProb;
Grid_new = zeros(size(Poly.total_pdf));
options = optimset('Display','off','TolX',StaticParams.criter_pdf);

wage = @(z_ag,kGrid) repmat((z_ag==0).*0.982*[2/3.06;4.12/3.06]...
                            +(z_ag==1).*1.054*[2/3.06;4.12/3.06],[1,length(kGrid)]);
wealth = @(z_ag,kGrid) wage(z_ag,kGrid)+repmat(kGrid,[2,1]);

Res.stErr_PPA1 = [];
Res.stErr_PPA2 = [];
Res.stErr_PPA3 = [];
Res.stErr_PPA4 = [];

Res.stErr_PFI1 = [];
Res.stErr_PFI2 = [];
Res.stErr_PFI3 = [];
Res.stErr_PFI4 = [];

Res.dynErr_PPA1 = [];
Res.dynErr_PPA2 = [];
Res.dynErr_PPA3 = [];
Res.dynErr_PPA4 = [];

Res.dynErr_PFI1 = [];
Res.dynErr_PFI2 = [];
Res.dynErr_PFI3 = [];
Res.dynErr_PFI4 = [];

Res.stErr_PPA1_algoGrid = [];
Res.stErr_PPA2_algoGrid = [];
Res.stErr_PPA3_algoGrid = [];
Res.stErr_PPA4_algoGrid = [];

Res.stErr_PFI1_algoGrid = [];
Res.stErr_PFI2_algoGrid = [];
Res.stErr_PFI3_algoGrid = [];
Res.stErr_PFI4_algoGrid = [];

Res.dynErr_PPA1_algoGrid = [];
Res.dynErr_PPA2_algoGrid = [];
Res.dynErr_PPA3_algoGrid = [];
Res.dynErr_PPA4_algoGrid = [];

Res.dynErr_PFI1_algoGrid = [];
Res.dynErr_PFI2_algoGrid = [];
Res.dynErr_PFI3_algoGrid = [];
Res.dynErr_PFI4_algoGrid = [];

Res.stErr_PPA1_pr = [];
Res.stErr_PPA2_pr = [];
Res.stErr_PPA3_pr = [];
Res.stErr_PPA4_pr = [];

Res.stErr_PFI1_pr = [];
Res.stErr_PFI2_pr = [];
Res.stErr_PFI3_pr = [];
Res.stErr_PFI4_pr = [];

Res.dynErr_PPA1_pr = [];
Res.dynErr_PPA2_pr = [];
Res.dynErr_PPA3_pr = [];
Res.dynErr_PPA4_pr = [];

Res.dynErr_PFI1_pr = [];
Res.dynErr_PFI2_pr = [];
Res.dynErr_PFI3_pr = [];
Res.dynErr_PFI4_pr = [];

%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------
for i=1:Res.nsamples
    % Euler Equ. Errors for PPA____________________________________________
    for method=0:1 
    for order_ct=1:M
        if method==0
            Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PPA.mat'));
            Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PPA.mat'));
        else
            Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
            Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'));
        end
        
        % Algorithm-implied policies
        K(method+1,order_ct,i) = StaticParams.kGrid*squeeze(pdf(method+1,order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
        sd(method+1,order_ct,i) = sqrt(StaticParams.kGrid.^2*squeeze(pdf(method+1,order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]);
        skew(method+1,order_ct,i) = ((StaticParams.kGrid-squeeze(K(method+1,order_ct,i))).^3 ...
                                      *squeeze(pdf(method+1,order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                                      ./squeeze(sd(method+1,order_ct,i))^3;
        kurt(method+1,order_ct,i) = ((StaticParams.kGrid-squeeze(K(method+1,order_ct,i))).^4 ...
                                      *squeeze(pdf(method+1,order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                                      ./squeeze(sd(method+1,order_ct,i))^4;
        k_prime = interp1(StaticParams.kGrid_pol,squeeze(...
                  [getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct,1:(Poly.phi_N+1)))')])',Res.kGrid)';
        c_prime = interp1(StaticParams.kGrid_pol,squeeze(...
                  [getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct,1:(Poly.phi_N+1)))')])',Res.kGrid)';   
        pr(method+1,order_ct,i) = getCApprox(Sol.distrGrid,squeeze(Sol.p(:,Res.z_ag(i)+1)),Sol.idx,squeeze(weights(method+1,order_ct,1:(Poly.phi_N+1))));
        if i==1
            kb_0 = kb(method+1,order_ct,1);
            kg_0 = kg(method+1,order_ct,1);
        else
            kb_0 = max(StaticParams.k_min,kb(method+1,order_ct,i-1));
            kg_0 = max(StaticParams.k_min,kg(method+1,order_ct,i-1));          
        end
        kb(method+1,order_ct,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct,1:(Poly.phi_N+1)))'))',kb_0)';
        kg(method+1,order_ct,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct,1:(Poly.phi_N+1)))'))',kg_0)';
        cb(method+1,order_ct,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(method+1,order_ct,1:(Poly.phi_N+1)))'))',kb_0)';   
        cg(method+1,order_ct,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(method+1,order_ct,1:(Poly.phi_N+1)))'))',kg_0)';   
        
        % standard Euler policies
        euler_weights = calcWeights(squeeze(cdf(method+1,order_ct,:,:)),Poly,StaticParams,Res.z_ag(i),...
                            interp1(Res.kGrid,sort(k_prime,2)',StaticParams.kGrid)');
        euler_expec = @(grid) (interp1(StaticParams.kGrid_pol,...
                               [getCApprox(Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')')...
                      .^(-StaticParams.gamma);
        euler_c_prime_aux = [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(k_prime(1,:)));...
                             (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(k_prime(2,:)))];
        cdf_aux = Poly.total_pdf;
        for b=1:size(cdf_aux,1)
            if b==1 || (b==2 && Res.z_ag(i)==0)
                cdf_aux(b,:,:) = cumsum(cdf_aux(b,:,:),3); 
                cdf_aux(b,2:end,:) = cdf_aux(b,2:end,:) + repmat(cumsum(cdf_aux(b,1:end-1,end)),[1,1,size(cdf_aux,3)]); 
            else
                cdf_aux(b,:,:) = cumsum(cdf_aux(b,:,:),2); 
                cdf_aux(b,:,2:end) = cdf_aux(b,:,2:end) + repmat(cumsum(cdf_aux(b,end,1:end-1)),[1,size(cdf_aux,3),1]); 
            end
            cdf_aux(b,:,:) = min(1,cdf_aux(b,:,:)./max(max(cdf_aux(b,:,:))));
        end
        Grid_new(1,:,:,:) = projectDistr(StaticParams.kGrid,squeeze(cdf(method+1,order_ct,1,:))',cdf_aux(1,:,:));
        Grid_new(2,:,:,:) = projectDistr(StaticParams.kGrid,squeeze(cdf(method+1,order_ct,Res.z_ag(i)+1,:))',cdf_aux(2,:,:));
        Grid_new(3,:,:,:) = projectDistr(StaticParams.kGrid,squeeze(cdf(method+1,order_ct,2,:))',cdf_aux(3,:,:));
        k_prime_aux = @(p) max(StaticParams.k_min,(wealth(Res.z_ag(i),Res.kGrid)-(euler_c_prime_aux./p).^(-1/StaticParams.gamma))./p);
        k_prime_aux2 = @(p) interpn(1:2,Res.kGrid,k_prime_aux(p),repmat((1:2)',[1,length(StaticParams.kGrid_pol)]),repmat(StaticParams.kGrid_pol,[2,1]));
        euler_p = ((wage(Res.z_ag(i),1)'*StaticParams.p(2*Res.z_ag(i)+(1:2))*2)/...
                   (sum(interp1(Res.kGrid,euler_c_prime_aux'.^(-1/StaticParams.gamma),StaticParams.kGrid)'.*squeeze(pdf(method+1,order_ct,:,:)),2)'*StaticParams.p(2*Res.z_ag(i)+(1:2))*2))...
                   ^StaticParams.gamma;
        euler_p = fzero(@(p) AggrAssetAlloc(p,k_prime_aux2(p),Grid_new,StaticParams,p,Res.z_ag(i),false,Poly),euler_p,options);
        euler_c_prime = wealth(Res.z_ag(i),Res.kGrid)-k_prime_aux(euler_p).*euler_p;
        
        % standard error
        Err = (c_prime-euler_c_prime)./max(eps,euler_c_prime);
        Err_algoGrid = (c_prime(:,algo_idx)-euler_c_prime(:,algo_idx))./max(eps,euler_c_prime(:,algo_idx));
        Err = Err(logical(squeeze(idx_posProb(method+1,order_ct,:,:))));
        Err_algoGrid = Err_algoGrid(logical(squeeze(idx_posProb(method+1,order_ct,:,algo_idx))));
        Err_pr = (pr(method+1,order_ct,i)-euler_p)./max(eps,euler_p);
        switch method
            case 0
                switch order_ct
                    case 1
                    Res.stErr_PPA1 = [Res.stErr_PPA1;Err];
                    Res.stErr_PPA1_algoGrid = [Res.stErr_PPA1_algoGrid;Err_algoGrid];
                    Res.stErr_PPA1_pr = [Res.stErr_PPA1_pr;Err_pr];
                    case 2
                    Res.stErr_PPA2 = [Res.stErr_PPA2;Err];
                    Res.stErr_PPA2_algoGrid = [Res.stErr_PPA2_algoGrid;Err_algoGrid];
                    Res.stErr_PPA2_pr = [Res.stErr_PPA2_pr;Err_pr];
                    case 3
                    Res.stErr_PPA3 = [Res.stErr_PPA3;Err];
                    Res.stErr_PPA3_algoGrid = [Res.stErr_PPA3_algoGrid;Err_algoGrid];
                    Res.stErr_PPA3_pr = [Res.stErr_PPA3_pr;Err_pr];
                    case 4
                    Res.stErr_PPA4 = [Res.stErr_PPA4;Err];
                    Res.stErr_PPA4_algoGrid = [Res.stErr_PPA4_algoGrid;Err_algoGrid];
                    Res.stErr_PPA4_pr = [Res.stErr_PPA4_pr;Err_pr];
                end
            case 1
                switch order_ct
                    case 1
                    Res.stErr_PFI1 = [Res.stErr_PFI1;Err];
                    Res.stErr_PFI1_algoGrid = [Res.stErr_PFI1_algoGrid;Err_algoGrid];
                    Res.stErr_PFI1_pr = [Res.stErr_PFI1_pr;Err_pr];
                    case 2
                    Res.stErr_PFI2 = [Res.stErr_PFI2;Err];
                    Res.stErr_PFI2_algoGrid = [Res.stErr_PFI2_algoGrid;Err_algoGrid];
                    Res.stErr_PFI2_pr = [Res.stErr_PFI2_pr;Err_pr];
                    case 3
                    Res.stErr_PFI3 = [Res.stErr_PFI3;Err];
                    Res.stErr_PFI3_algoGrid = [Res.stErr_PFI3_algoGrid;Err_algoGrid];
                    Res.stErr_PFI3_pr = [Res.stErr_PFI3_pr;Err_pr];
                    case 4
                    Res.stErr_PFI4 = [Res.stErr_PFI4;Err];
                    Res.stErr_PFI4_algoGrid = [Res.stErr_PFI4_algoGrid;Err_algoGrid];
                    Res.stErr_PFI4_pr = [Res.stErr_PFI4_pr;Err_pr];
                end
        end
        % update distribution
        if i<Res.nsamples
            pdf(method+1,order_ct,:,:) = calcEoPDistr(interp1(Res.kGrid,sort(k_prime,2)',StaticParams.kGrid)',...
                           squeeze(cdf(method+1,order_ct,:,:)),StaticParams.kGrid);
            pdf(method+1,order_ct,:,:) = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                           .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                           *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(pdf,4)])...
                           .*squeeze(pdf(method+1,order_ct,:,:)));
            idx = max(find(squeeze(pdf(method+1,order_ct,1,:))>eps,1,'last'),find(squeeze(pdf(method+1,order_ct,2,:))>eps,1,'last'));
            idx_posProb(method+1,order_ct,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
            pdf = pdf.*(pdf>eps);
            pdf = pdf./repmat(sum(pdf,4),[1,1,1,size(pdf,4)]);
            cdf = cummax(min(1,cumsum(pdf,4)),4);
            cdf(method+1,order_ct,1,cdf(method+1,order_ct,1,:)==max(cdf(method+1,order_ct,1,:))) = 1;
            cdf(method+1,order_ct,2,cdf(method+1,order_ct,2,:)==max(cdf(method+1,order_ct,2,:))) = 1;
            weights(method+1,order_ct,1:(Poly.phi_N+1))...
            = calcWeights(squeeze(cdf(method+1,order_ct,:,:)),Poly,StaticParams,Res.z_ag(i+1));
        end
        
        % dynamic Euler policies and distributions
        dyn_euler_K(method+1,order_ct,i) = StaticParams.kGrid*squeeze(dyn_pdf(method+1,order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
        k_prime_next = [getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(dyn_weights(method+1,order_ct,1:(Poly.phi_N+1)))');...
                        getCApprox(Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(dyn_weights(method+1,order_ct,1:(Poly.phi_N+1)))')];  
        euler_weights = calcWeights(squeeze(dyn_cdf(method+1,order_ct,:,:)),Poly,StaticParams,Res.z_ag(i),...
                        interp1(StaticParams.kGrid_pol,sort(k_prime_next,2)',StaticParams.kGrid)');
        euler_expec = @(grid) (interp1(StaticParams.kGrid_pol,...
                               [getCApprox(Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')')...
                      .^(-StaticParams.gamma);
        euler_c_prime_aux = [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*euler_expec(interp1(StaticParams.kGrid_pol,sort(k_prime_next(1,:),2)',Res.kGrid)'));...
                             (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*euler_expec(interp1(StaticParams.kGrid_pol,sort(k_prime_next(2,:),2)',Res.kGrid)'))];
        Grid_new(1,:,:,:) = projectDistr(StaticParams.kGrid,squeeze(cdf(method+1,order_ct,1,:))',cdf_aux(1,:,:));
        Grid_new(2,:,:,:) = projectDistr(StaticParams.kGrid,squeeze(cdf(method+1,order_ct,Res.z_ag(i)+1,:))',cdf_aux(2,:,:));
        Grid_new(3,:,:,:) = projectDistr(StaticParams.kGrid,squeeze(cdf(method+1,order_ct,2,:))',cdf_aux(3,:,:));
        dyn_euler_k_prime = @(p) max(StaticParams.k_min,(wealth(Res.z_ag(i),Res.kGrid)-(euler_c_prime_aux./p).^(-1/StaticParams.gamma))./p);
        k_prime_aux2 = @(p) interpn(1:2,Res.kGrid,dyn_euler_k_prime(p),repmat((1:2)',[1,length(StaticParams.kGrid_pol)]),repmat(StaticParams.kGrid_pol,[2,1]));
        dyn_euler_p = ((wage(Res.z_ag(i),1)'*StaticParams.p(2*Res.z_ag(i)+(1:2))*2)/...
                   (sum(interp1(Res.kGrid,euler_c_prime_aux'.^(-1/StaticParams.gamma),StaticParams.kGrid)'.*squeeze(dyn_pdf(method+1,order_ct,:,:)),2)'*StaticParams.p(2*Res.z_ag(i)+(1:2))*2))...
                   ^StaticParams.gamma;
        dyn_euler_p = fzero(@(p) AggrAssetAlloc(p,k_prime_aux2(p),Grid_new,StaticParams,p,Res.z_ag(i),false,Poly),dyn_euler_p,options);
        dyn_euler_k_prime = dyn_euler_k_prime(dyn_euler_p);
        dyn_euler_c_prime = wealth(Res.z_ag(i),Res.kGrid)-dyn_euler_k_prime.*dyn_euler_p;
        
        % dynamic error
        Err = (c_prime-dyn_euler_c_prime)./max(eps,dyn_euler_c_prime);
        Err_algoGrid = (c_prime(:,algo_idx)-dyn_euler_c_prime(:,algo_idx))./max(eps,dyn_euler_c_prime(:,algo_idx));
        Err = Err(logical(squeeze(dyn_idx_posProb(method+1,order_ct,:,:))));
        Err_algoGrid = Err_algoGrid(logical(squeeze(dyn_idx_posProb(method+1,order_ct,:,algo_idx))));
        Err_pr = (pr(method+1,order_ct,i)-dyn_euler_p)./max(eps,dyn_euler_p);
        switch method
            case 0
                switch order_ct
                    case 1
                    Res.dynErr_PPA1 = [Res.dynErr_PPA1;Err];
                    Res.dynErr_PPA1_algoGrid = [Res.dynErr_PPA1_algoGrid;Err_algoGrid];
                    Res.dynErr_PPA1_pr = [Res.dynErr_PPA1_pr;Err_pr];
                    case 2
                    Res.dynErr_PPA2 = [Res.dynErr_PPA2;Err];
                    Res.dynErr_PPA2_algoGrid = [Res.dynErr_PPA2_algoGrid;Err_algoGrid];
                    Res.dynErr_PPA2_pr = [Res.dynErr_PPA2_pr;Err_pr];
                    case 3
                    Res.dynErr_PPA3 = [Res.dynErr_PPA3;Err];
                    Res.dynErr_PPA3_algoGrid = [Res.dynErr_PPA3_algoGrid;Err_algoGrid];
                    Res.dynErr_PPA3_pr = [Res.dynErr_PPA3_pr;Err_pr];
                    case 4
                    Res.dynErr_PPA4 = [Res.dynErr_PPA4;Err];
                    Res.dynErr_PPA4_algoGrid = [Res.dynErr_PPA4_algoGrid;Err_algoGrid];
                    Res.dynErr_PPA4_pr = [Res.dynErr_PPA4_pr;Err_pr];
                end
            case 1
                switch order_ct
                    case 1
                    Res.dynErr_PFI1 = [Res.dynErr_PFI1;Err];
                    Res.dynErr_PFI1_algoGrid = [Res.dynErr_PFI1_algoGrid;Err_algoGrid];
                    Res.dynErr_PFI1_pr = [Res.dynErr_PFI1_pr;Err_pr];
                    case 2
                    Res.dynErr_PFI2 = [Res.dynErr_PFI2;Err];
                    Res.dynErr_PFI2_algoGrid = [Res.dynErr_PFI2_algoGrid;Err_algoGrid];
                    Res.dynErr_PFI2_pr = [Res.dynErr_PFI2_pr;Err_pr];
                    case 3
                    Res.dynErr_PFI3 = [Res.dynErr_PFI3;Err];
                    Res.dynErr_PFI3_algoGrid = [Res.dynErr_PFI3_algoGrid;Err_algoGrid];
                    Res.dynErr_PFI3_pr = [Res.dynErr_PFI3_pr;Err_pr];
                    case 4
                    Res.dynErr_PFI4 = [Res.dynErr_PFI4;Err];
                    Res.dynErr_PFI4_algoGrid = [Res.dynErr_PFI4_algoGrid;Err_algoGrid];
                    Res.dynErr_PFI4_pr = [Res.dynErr_PFI4_pr;Err_pr];
                end
        end
        % update distribution
        if i<Res.nsamples
            dyn_pdf(method+1,order_ct,:,:) = calcEoPDistr(interp1(Res.kGrid,sort(dyn_euler_k_prime,2)',StaticParams.kGrid)',...
                           squeeze(dyn_cdf(method+1,order_ct,:,:)),StaticParams.kGrid);
            dyn_pdf(method+1,order_ct,:,:) = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                           .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                           *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(dyn_pdf,4)])...
                           .*squeeze(dyn_pdf(method+1,order_ct,:,:)));
            idx = max(find(squeeze(dyn_pdf(method+1,order_ct,1,:))>eps,1,'last'),find(squeeze(dyn_pdf(method+1,order_ct,2,:))>eps,1,'last'));
            dyn_idx_posProb(method+1,order_ct,:,(Res.kGrid<=StaticParams.kGrid(idx))) = 1;
            dyn_pdf = dyn_pdf.*(dyn_pdf>eps);
            dyn_pdf = dyn_pdf./repmat(sum(dyn_pdf,4),[1,1,1,size(dyn_pdf,4)]);
            dyn_cdf = cummax(min(1,cumsum(dyn_pdf,4)),4);
            dyn_cdf(method+1,order_ct,1,dyn_cdf(method+1,order_ct,1,:)==max(dyn_cdf(method+1,order_ct,1,:))) = 1;
            dyn_cdf(method+1,order_ct,2,dyn_cdf(method+1,order_ct,2,:)==max(dyn_cdf(method+1,order_ct,2,:))) = 1;
            dyn_weights(method+1,order_ct,1:(Poly.phi_N+1))...
            = calcWeights(squeeze(dyn_cdf(method+1,order_ct,:,:)),Poly,StaticParams,Res.z_ag(i+1));
        end
    end
    end
    
    if i<Res.nsamples    
        if Res.z_ag(i)==0
            fill(x,y,'r','EdgeColor','none')
        else
            fill(x,y,'g','EdgeColor','none')
        end
        set(gca,'nextplot','add');
        pdf_aux = [squeeze(pdf(1,M,1,:))',squeeze(pdf(1,M,1,2:end))'];
        pdf_aux = pdf_aux(idx_pdf);
        plot(grid,pdf_aux,'m','LineWidth',2);
        set(gca,'nextplot','add');
        pdf_aux = [squeeze(pdf(1,M,2,:))',squeeze(pdf(1,M,2,2:end))'];
        pdf_aux = pdf_aux(idx_pdf);
        plot(grid,pdf_aux,'b','LineWidth',2);
        set(gca,'nextplot','replacechildren');
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
close(v)

Res.kb = kb;
Res.kg = kg;
Res.cb = cb;
Res.cg = cg;

Res.pr = pr;
Res.pr_mean = mean(Res.pr,3);
Res.pr_mean0 = mean(Res.pr(:,:,(Res.z_ag==0)),3);
Res.pr_mean1 = mean(Res.pr(:,:,(Res.z_ag==1)),3);
Res.pr_sd = std(Res.pr,[],3);
Res.pr_sd0 = std(Res.pr(:,:,(Res.z_ag==0)),[],3);
Res.pr_sd1 = std(Res.pr(:,:,(Res.z_ag==1)),[],3);
Res.pr_skew = skewness(Res.pr,[],3);
Res.pr_skew0 = skewness(Res.pr(:,:,(Res.z_ag==0)),[],3);
Res.pr_skew1 = skewness(Res.pr(:,:,(Res.z_ag==1)),[],3);
Res.pr_kurt = kurtosis(Res.pr,[],3);
Res.pr_kurt0 = kurtosis(Res.pr(:,:,(Res.z_ag==0)),[],3);
Res.pr_kurt1 = kurtosis(Res.pr(:,:,(Res.z_ag==1)),[],3);

Res.K = K;
Res.sd = sd;
Res.skew = skew;
Res.kurt = kurt;
Res.dyn_euler_K = dyn_euler_K;

%--------------------------------------------------------------------------
% Compute the Euler equation error statistics
%--------------------------------------------------------------------------
Res.median_stErr_PPA = [median(abs(Res.stErr_PPA1));median(abs(Res.stErr_PPA2));median(abs(Res.stErr_PPA3));median(abs(Res.stErr_PPA4))];
Res.median_stErr_PFI = [median(abs(Res.stErr_PFI1));median(abs(Res.stErr_PFI2));median(abs(Res.stErr_PFI3));median(abs(Res.stErr_PFI4))];

Res.mean_stErr_PPA = [mean(abs(Res.stErr_PPA1));mean(abs(Res.stErr_PPA2));mean(abs(Res.stErr_PPA3));mean(abs(Res.stErr_PPA4))];
Res.mean_stErr_PFI = [mean(abs(Res.stErr_PFI1));mean(abs(Res.stErr_PFI2));mean(abs(Res.stErr_PFI3));mean(abs(Res.stErr_PFI4))];

Res.min_stErr_PPA = [min(abs(Res.stErr_PPA1));min(abs(Res.stErr_PPA2));min(abs(Res.stErr_PPA3));min(abs(Res.stErr_PPA4))];
Res.min_stErr_PFI = [min(abs(Res.stErr_PFI1));min(abs(Res.stErr_PFI2));min(abs(Res.stErr_PFI3));min(abs(Res.stErr_PFI4))];

Res.max_stErr_PPA = [max(abs(Res.stErr_PPA1));max(abs(Res.stErr_PPA2));max(abs(Res.stErr_PPA3));max(abs(Res.stErr_PPA4))];
Res.max_stErr_PFI = [max(abs(Res.stErr_PFI1));max(abs(Res.stErr_PFI2));max(abs(Res.stErr_PFI3));max(abs(Res.stErr_PFI4))];

Res.median_dynErr_PPA = [median(abs(Res.dynErr_PPA1));median(abs(Res.dynErr_PPA2));median(abs(Res.dynErr_PPA3));median(abs(Res.dynErr_PPA4))];
Res.median_dynErr_PFI = [median(abs(Res.dynErr_PFI1));median(abs(Res.dynErr_PFI2));median(abs(Res.dynErr_PFI3));median(abs(Res.dynErr_PFI4))];

Res.mean_dynErr_PPA = [mean(abs(Res.dynErr_PPA1));mean(abs(Res.dynErr_PPA2));mean(abs(Res.dynErr_PPA3));mean(abs(Res.dynErr_PPA4))];
Res.mean_dynErr_PFI = [mean(abs(Res.dynErr_PFI1));mean(abs(Res.dynErr_PFI2));mean(abs(Res.dynErr_PFI3));mean(abs(Res.dynErr_PFI4))];

Res.min_dynErr_PPA = [min(abs(Res.dynErr_PPA1));min(abs(Res.dynErr_PPA2));min(abs(Res.dynErr_PPA3));min(abs(Res.dynErr_PPA4))];
Res.min_dynErr_PFI = [min(abs(Res.dynErr_PFI1));min(abs(Res.dynErr_PFI2));min(abs(Res.dynErr_PFI3));min(abs(Res.dynErr_PFI4))];

Res.max_dynErr_PPA = [max(abs(Res.dynErr_PPA1));max(abs(Res.dynErr_PPA2));max(abs(Res.dynErr_PPA3));max(abs(Res.dynErr_PPA4))];
Res.max_dynErr_PFI = [max(abs(Res.dynErr_PFI1));max(abs(Res.dynErr_PFI2));max(abs(Res.dynErr_PFI3));max(abs(Res.dynErr_PFI4))];


Res.median_stErr_PPA_pr = [median(abs(Res.stErr_PPA1_pr));median(abs(Res.stErr_PPA2_pr));median(abs(Res.stErr_PPA3_pr));median(abs(Res.stErr_PPA4_pr))];
Res.median_stErr_PFI_pr = [median(abs(Res.stErr_PFI1_pr));median(abs(Res.stErr_PFI2_pr));median(abs(Res.stErr_PFI3_pr));median(abs(Res.stErr_PFI4_pr))];

Res.mean_stErr_PPA_pr = [mean(abs(Res.stErr_PPA1_pr));mean(abs(Res.stErr_PPA2_pr));mean(abs(Res.stErr_PPA3_pr));mean(abs(Res.stErr_PPA4_pr))];
Res.mean_stErr_PFI_pr = [mean(abs(Res.stErr_PFI1_pr));mean(abs(Res.stErr_PFI2_pr));mean(abs(Res.stErr_PFI3_pr));mean(abs(Res.stErr_PFI4_pr))];

Res.min_stErr_PPA_pr = [min(abs(Res.stErr_PPA1_pr));min(abs(Res.stErr_PPA2_pr));min(abs(Res.stErr_PPA3_pr));min(abs(Res.stErr_PPA4_pr))];
Res.min_stErr_PFI_pr = [min(abs(Res.stErr_PFI1_pr));min(abs(Res.stErr_PFI2_pr));min(abs(Res.stErr_PFI3_pr));min(abs(Res.stErr_PFI4_pr))];

Res.max_stErr_PPA_pr = [max(abs(Res.stErr_PPA1_pr));max(abs(Res.stErr_PPA2_pr));max(abs(Res.stErr_PPA3_pr));max(abs(Res.stErr_PPA4_pr))];
Res.max_stErr_PFI_pr = [max(abs(Res.stErr_PFI1_pr));max(abs(Res.stErr_PFI2_pr));max(abs(Res.stErr_PFI3_pr));max(abs(Res.stErr_PFI4_pr))];

Res.median_dynErr_PPA_pr = [median(abs(Res.dynErr_PPA1_pr));median(abs(Res.dynErr_PPA2_pr));median(abs(Res.dynErr_PPA3_pr));median(abs(Res.dynErr_PPA4_pr))];
Res.median_dynErr_PFI_pr = [median(abs(Res.dynErr_PFI1_pr));median(abs(Res.dynErr_PFI2_pr));median(abs(Res.dynErr_PFI3_pr));median(abs(Res.dynErr_PFI4_pr))];

Res.mean_dynErr_PPA_pr = [mean(abs(Res.dynErr_PPA1_pr));mean(abs(Res.dynErr_PPA2_pr));mean(abs(Res.dynErr_PPA3_pr));mean(abs(Res.dynErr_PPA4_pr))];
Res.mean_dynErr_PFI_pr = [mean(abs(Res.dynErr_PFI1_pr));mean(abs(Res.dynErr_PFI2_pr));mean(abs(Res.dynErr_PFI3_pr));mean(abs(Res.dynErr_PFI4_pr))];

Res.min_dynErr_PPA_pr = [min(abs(Res.dynErr_PPA1_pr));min(abs(Res.dynErr_PPA2_pr));min(abs(Res.dynErr_PPA3_pr));min(abs(Res.dynErr_PPA4_pr))];
Res.min_dynErr_PFI_pr = [min(abs(Res.dynErr_PFI1_pr));min(abs(Res.dynErr_PFI2_pr));min(abs(Res.dynErr_PFI3_pr));min(abs(Res.dynErr_PFI4_pr))];

Res.max_dynErr_PPA_pr = [max(abs(Res.dynErr_PPA1_pr));max(abs(Res.dynErr_PPA2_pr));max(abs(Res.dynErr_PPA3_pr));max(abs(Res.dynErr_PPA4_pr))];
Res.max_dynErr_PFI_pr = [max(abs(Res.dynErr_PFI1_pr));max(abs(Res.dynErr_PFI2_pr));max(abs(Res.dynErr_PFI3_pr));max(abs(Res.dynErr_PFI4_pr))];

save(strcat(folder,'/Res.mat'),'-struct','Res');
end