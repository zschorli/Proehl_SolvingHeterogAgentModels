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
function results_Errors(StaticParams,folder)

% truncation order 
M = 3;

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
preSol = load(strcat(folder,'/Sol',num2str(3),'_PFI.mat'));
pdf_start = preSol.pdf_cond(1:2,:);
cdf_start = cummax(min(1,cumsum(pdf_start,2)),2);
cdf_start(1,cdf_start(1,:)==max(cdf_start(1,:)))=1;
cdf_start(2,cdf_start(2,:)==max(cdf_start(2,:)))=1;
pdf_start_l = preSol.pdf_cond_l(1:2,:);

% Set variables for PFI order 0 to M
pdf = repmat(permute(pdf_start,[3,1,2]),[M+1,1,1]);
cdf = repmat(permute(cdf_start,[3,1,2]),[M+1,1,1]);
weights = zeros(M+1,M+1);
for order_ct=0:M
    Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
    weights(order_ct+1,1:(Poly.phi_N+1))...
    = calcWeights(squeeze(cdf(order_ct+1,:,:)),Poly,StaticParams,0);
end
pdf_l = repmat(permute(pdf_start_l,[3,1,2]),[M+1,1,1]);
kb = zeros(M+1,Res.nsamples);
kg = zeros(M+1,Res.nsamples);
cb = zeros(M+1,Res.nsamples);
cg = zeros(M+1,Res.nsamples);
lb = zeros(M+1,Res.nsamples);
lg = zeros(M+1,Res.nsamples);
K = zeros(M+1,Res.nsamples);
L = zeros(M+1,Res.nsamples);
sd = zeros(M+1,Res.nsamples);
skew = zeros(M+1,Res.nsamples);
kurt = zeros(M+1,Res.nsamples);
gini = zeros(M+1,Res.nsamples);
K_pr_LoM = zeros(M+1,Res.nsamples);
L_pr_LoM = zeros(M+1,Res.nsamples);

% Simulation video
v = VideoWriter(strcat(folder,'/Distr_Sim.avi'));
open(v);
k_max=30;
x=[0 0 k_max k_max];
y=[-5e-3 0 0 -5e-3];


%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------
for i=1:Res.nsamples
    % Euler Equ. Errors for PPA____________________________________________ 
    for order_ct=0:M
        Poly = load(strcat(folder,'/Poly',num2str(order_ct),'_PFI.mat'));
        Sol = load(strcat(folder,'/Sol',num2str(order_ct),'_PFI.mat'));

        % update labor supply distribution
        l_prime = interp1(StaticParams.kGrid_pol,squeeze(...
                  [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))')])',StaticParams.kGrid)';
        pdf_l(order_ct+1,:,:) = calcEoPDistr(1-l_prime,squeeze(cdf(order_ct+1,:,:)),StaticParams.lGrid);
        pdf_l(order_ct+1,:,:) = pdf_l(order_ct+1,:,size(pdf_l,3):-1:1);
        pdf_l(order_ct+1,:,:) = pdf_l(order_ct+1,:,:).*(pdf_l(order_ct+1,:,:)>eps);
        pdf_l(order_ct+1,:,:) = pdf_l(order_ct+1,:,:)./repmat(sum(pdf_l(order_ct+1,:,:),3),[1,1,size(pdf_l,3)]);

        % update video
        if order_ct == M
        figure(1)%case_nr
        subplot(1,2,1)
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
        pdf_plot = 2*StaticParams.p(Res.z_ag(i)*2+(1:2))'*squeeze(pdf(end,:,:));
        p3=plot(StaticParams.kGrid,pdf_plot,'b','LineWidth',2);
        set(gca,'nextplot','replacechildren');
        legend([p3,p4,p5],'capital distribution', 'boom', 'recession','Location','NorthEast','AutoUpdate','off')
        legend('boxoff')
        axis([0 k_max -8e-5 5e-3])
        xlabel('k')
        ylabel('probability')
        subplot(1,2,2)
        if Res.z_ag(i)==1
            p8=fill(x,y,'r','EdgeColor','none');
            set(gca,'nextplot','add');
            p7=fill(x,y,'g','EdgeColor','none');
        else
            p7=fill(x,y,'g','EdgeColor','none');
            set(gca,'nextplot','add');
            p8=fill(x,y,'r','EdgeColor','none');
        end
        set(gca,'nextplot','add');
        p6=plot(StaticParams.lGrid*squeeze(pdf_l(end,2,:)).*ones(1,2),[0 1],'m','LineWidth',2);
        set(gca,'nextplot','add');
        p9=plot(StaticParams.lGrid,squeeze(pdf_l(end,2,:)),'b','LineWidth',2);
        set(gca,'nextplot','replacechildren');
        legend([p6,p9,p7,p8],'aggr. labor supply','labor supply distr.', 'boom', 'recession','Location','NorthEast','AutoUpdate','off')
        legend('boxoff')
        axis([0.2 0.45 -1.5e-3 0.08])
        xlabel('l')
        ylabel('probability')
        frame = getframe(gcf);
        writeVideo(v,frame);
        end
        
        % Algorithm-implied policies
        L(order_ct+1,i) = StaticParams.lGrid*squeeze(pdf_l(order_ct+1,2,:));
        K(order_ct+1,i) = StaticParams.kGrid*squeeze(pdf(order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
        sd(order_ct+1,i) = sqrt(StaticParams.kGrid.^2*squeeze(pdf(order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]);
        skew(order_ct+1,i) = ((StaticParams.kGrid-squeeze(K(order_ct+1,i))).^3 ...
                                      *squeeze(pdf(order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                                      ./squeeze(sd(order_ct+1,i))^3;
        kurt(order_ct+1,i) = ((StaticParams.kGrid-squeeze(K(order_ct+1,i))).^4 ...
                                      *squeeze(pdf(order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])...
                                      ./squeeze(sd(order_ct+1,i))^4;
        gini(order_ct+1,i) = sum(sum((repmat(StaticParams.kGrid',[1,length(StaticParams.kGrid)])-repmat(StaticParams.kGrid,[length(StaticParams.kGrid),1]))...
                                      .*repmat((squeeze(pdf(order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)]),[1,length(StaticParams.kGrid)])...
                                      .*repmat((squeeze(pdf(order_ct+1,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])',[length(StaticParams.kGrid),1]),2),1)...
                                      ./(2*squeeze(K(order_ct+1,i)));
        k_prime = interp1(StaticParams.kGrid_pol,squeeze(...
                  [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))')])',Res.kGrid)';
        if i==1
            kb_0 = K(order_ct+1,i);
            kg_0 = K(order_ct+1,i);
        else
            kb_0 = max(StaticParams.k_min,kb(order_ct+1,i-1));
            kg_0 = max(StaticParams.k_min,kg(order_ct+1,i-1));        
        end
        kb(order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))'))',kb_0)';
        kg(order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.k_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))'))',kg_0)';
        cb(order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))'))',kb_0)';   
        cg(order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))'))',kg_0)';   
        lb(order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))'))',kb_0)';   
        lg(order_ct+1,i) = interp1(StaticParams.kGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))'))',kg_0)';   
        if i<Res.nsamples
            LoM = calcAggregates(interp1(Res.kGrid,k_prime',StaticParams.kGrid)',...
                  squeeze(sum(repmat(squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))',[1,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:))]).*Poly.total,1)),...
                  StaticParams,Poly);
            K_pr_LoM(order_ct+1,i) = LoM(Res.z_ag(i)+1);
            [~,~,~,LoM_L] = calcAggregates(l_prime,...
                  squeeze(sum(repmat(squeeze(weights(order_ct+1,1:(Poly.phi_N+1)))',[1,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:))]).*Poly.total,1)),...
                  StaticParams,Poly);
            L_pr_LoM(order_ct+1,i) = LoM_L(Res.z_ag(i)+1);
        end
        
        % update distribution
        if i<Res.nsamples
            pdf(order_ct+1,:,:) = calcEoPDistr(interp1(Res.kGrid,k_prime',StaticParams.kGrid)',...
                           squeeze(cdf(order_ct+1,:,:)),StaticParams.kGrid);
            pdf(order_ct+1,:,:) = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                           .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                           *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(pdf,3)])...
                           .*squeeze(pdf(order_ct+1,:,:)));
            pdf = pdf.*(pdf>eps);
            pdf = pdf./repmat(sum(pdf,3),[1,1,size(pdf,3)]);
            cdf = cummax(min(1,cumsum(pdf,3)),3);
            cdf(order_ct+1,1,cdf(order_ct+1,1,:)==max(cdf(order_ct+1,1,:))) = 1;
            cdf(order_ct+1,2,cdf(order_ct+1,2,:)==max(cdf(order_ct+1,2,:))) = 1;
            weights(order_ct+1,1:(Poly.phi_N+1))...
            = calcWeights(squeeze(cdf(order_ct+1,:,:)),Poly,StaticParams,Res.z_ag(i+1));
        end
    end
end
close(v)

Res.kb = kb;
Res.kg = kg;
Res.cb = cb;
Res.cg = cg;

Res.K = K;
Res.L = L;
Res.sd = sd;
Res.skew = skew;
Res.kurt = kurt;
Res.gini = gini;
Res.K_pr_LoM = K_pr_LoM;
Res.L_pr_LoM = L_pr_LoM;


Res.LoMErr_PFI0 = squeeze(log10(abs(Res.K_pr_LoM(1,1:end-1)./Res.K(1,2:end)-1)));
Res.LoMErr_PFI1 = squeeze(log10(abs(Res.K_pr_LoM(2,1:end-1)./Res.K(2,2:end)-1)));
Res.LoMErr_PFI2 = squeeze(log10(abs(Res.K_pr_LoM(3,1:end-1)./Res.K(3,2:end)-1)));
Res.LoMErr_PFI3 = squeeze(log10(abs(Res.K_pr_LoM(4,1:end-1)./Res.K(4,2:end)-1)));
if M>3
    Res.LoMErr_PFI4 = squeeze(log10(abs(Res.K_pr_LoM(5,1:end-1)./Res.K(5,2:end)-1)));
end
Res.LoMErrL_PFI0 = squeeze(log10(abs(Res.L_pr_LoM(1,1:end-1)./Res.L(1,2:end)-1)));
Res.LoMErrL_PFI1 = squeeze(log10(abs(Res.L_pr_LoM(2,1:end-1)./Res.L(2,2:end)-1)));
Res.LoMErrL_PFI2 = squeeze(log10(abs(Res.L_pr_LoM(3,1:end-1)./Res.L(3,2:end)-1)));
Res.LoMErrL_PFI3 = squeeze(log10(abs(Res.L_pr_LoM(4,1:end-1)./Res.L(4,2:end)-1)));
if M>3
    Res.LoMErrL_PFI4 = squeeze(log10(abs(Res.L_pr_LoM(5,1:end-1)./Res.L(5,2:end)-1)));
end
Res.LoMErrKL_PFI0 = squeeze(log10(abs((Res.K_pr_LoM(1,1:end-1)./Res.L_pr_LoM(1,1:end-1))./(Res.K(1,2:end)./Res.L(1,2:end))-1)));
Res.LoMErrKL_PFI1 = squeeze(log10(abs((Res.K_pr_LoM(2,1:end-1)./Res.L_pr_LoM(2,1:end-1))./(Res.K(2,2:end)./Res.L(2,2:end))-1)));
Res.LoMErrKL_PFI2 = squeeze(log10(abs((Res.K_pr_LoM(3,1:end-1)./Res.L_pr_LoM(3,1:end-1))./(Res.K(3,2:end)./Res.L(3,2:end))-1)));
Res.LoMErrKL_PFI3 = squeeze(log10(abs((Res.K_pr_LoM(4,1:end-1)./Res.L_pr_LoM(4,1:end-1))./(Res.K(4,2:end)./Res.L(4,2:end))-1)));
if M>3
    Res.LoMErrKL_PFI4 = squeeze(log10(abs((Res.K_pr_LoM(5,1:end-1)./Res.L_pr_LoM(5,1:end-1))./(Res.K(5,2:end)./Res.L(5,2:end))-1)));
end

save(strcat(folder,'/Res.mat'),'-struct','Res');
end