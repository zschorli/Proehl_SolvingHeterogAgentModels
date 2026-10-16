% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This file computes the expected ergodic distributions of the different 
% algorithm solutions.
%__________________________________________________________________________
function results_Errors(StaticParams,folder,order)

% truncation order 
M = order+1;

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
Res.aGrid = unique([StaticParams.aGrid_pol,linspace(StaticParams.a_min,StaticParams.a_max,1001)]);
algo_idx = logical(sum(repmat(Res.aGrid,[length(StaticParams.aGrid_pol),1])==repmat(StaticParams.aGrid_pol',[1,length(Res.aGrid)]),1));
Sol = load(strcat(folder,'/Sol3_PFI.mat'));
pdf_start = 0.5.*(Sol.pdf_cond(1:2,:)+Sol.pdf_cond(3:4,:));
cdf_start = Sol.cdf;

v = VideoWriter(strcat(folder,'/Distr_Sim.avi'));
open(v);
a_max=0.3*StaticParams.a_max;
x=[StaticParams.a_min StaticParams.a_min a_max a_max];
y=[-0.015 0 0 -0.015];
[grid,idx_pdf] = sort([StaticParams.aGrid,StaticParams.aGrid(1:end-1)+1e-5]);
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
axis([StaticParams.a_min a_max -0.0005 0.02])
title('Simulation of Net Worth Distributions')
xlabel('a')
ylabel('probability')
frame = getframe(gcf);
writeVideo(v,frame);

% Set variables for PCE order 1 to M
pdf = repmat(permute(pdf_start,[3,1,2]),[M,1,1]);
cdf = repmat(permute(cdf_start,[3,1,2]),[M,1,1]);
weights = zeros(M,M);
for order_ct=1:M
    Poly = load(strcat(folder,'/Poly',num2str(order_ct-1),'_PFI.mat'));
    weights(order_ct,1:(Poly.phi_N+1))...
    = calcWeights(squeeze(cdf(order_ct,:,:)),Poly,StaticParams,0);
end
dyn_pdf = pdf;
dyn_cdf = cdf;
dyn_weights = weights;
ab = zeros(M,Res.nsamples,1);
ag = zeros(M,Res.nsamples,1);
hb = zeros(M,Res.nsamples,1);
hg = zeros(M,Res.nsamples,1);
bb = zeros(M,Res.nsamples,1);
bg = zeros(M,Res.nsamples,1);
cb = zeros(M,Res.nsamples,1);
cg = zeros(M,Res.nsamples,1);
q = zeros(M,Res.nsamples,1);
r = zeros(M,Res.nsamples,1);
A = zeros(M,Res.nsamples,1);
sd = zeros(M,Res.nsamples,1);
skew = zeros(M,Res.nsamples,1);
kurt = zeros(M,Res.nsamples,1);
gini = zeros(M,Res.nsamples,1);
H = zeros(M,Res.nsamples,1);
Borrow = zeros(M,Res.nsamples,1);
Lend = zeros(M,Res.nsamples,1);
dyn_euler_A = zeros(M,Res.nsamples,1);
dyn_euler_H = zeros(M,Res.nsamples,1);
dyn_euler_Borrow = zeros(M,Res.nsamples,1);
dyn_euler_Lend = zeros(M,Res.nsamples,1);
idx_posProb = zeros(M,2,length(Res.aGrid));
idx = max(find(pdf(1,1,:)>eps,1,'last'),find(pdf(1,2,:)>eps,1,'last'));
idx_posProb(:,:,(Res.aGrid<=StaticParams.aGrid(idx))) = 1;
dyn_idx_posProb = idx_posProb;

innerAggr = @(c,h) c.^((StaticParams.eps-1)/StaticParams.eps)...
               +StaticParams.eta^(1/StaticParams.eps).*...
               h.^((StaticParams.eps-1)/StaticParams.eps);
Aggr = @(c,h) (innerAggr(c,h)).^(StaticParams.eps/(StaticParams.eps-1));
Aggr_derC = @(c,h) (innerAggr(c,h)).^(StaticParams.eps/(StaticParams.eps-1)-1)...
               .*c.^(-1/StaticParams.eps);
F_uq = @(r,q,q_pr) max(eps,(1+r).*q-q_pr);
a_large = repmat(Res.aGrid,[2,1]);
h_bound = @(q) 1/(1-StaticParams.theta).*a_large./q;
F_h = @(r,q,q_pr,c) min((StaticParams.eta.*(F_uq(r,q,q_pr)).^(-StaticParams.eps).*c)...
      ,h_bound(q));

Res.stErr_PFI1 = [];
Res.stErr_PFI2 = [];
Res.stErr_PFI3 = [];
Res.stErr_PFI4 = [];

Res.dynErr_PFI1 = [];
Res.dynErr_PFI2 = [];
Res.dynErr_PFI3 = [];
Res.dynErr_PFI4 = [];

Res.stErr_PFI1_algoGrid = [];
Res.stErr_PFI2_algoGrid = [];
Res.stErr_PFI3_algoGrid = [];
Res.stErr_PFI4_algoGrid = [];

Res.dynErr_PFI1_algoGrid = [];
Res.dynErr_PFI2_algoGrid = [];
Res.dynErr_PFI3_algoGrid = [];
Res.dynErr_PFI4_algoGrid = [];

Res.stErr_PFI1_q = [];
Res.stErr_PFI2_q = [];
Res.stErr_PFI3_q = [];
Res.stErr_PFI4_q = [];

Res.dynErr_PFI1_q = [];
Res.dynErr_PFI2_q = [];
Res.dynErr_PFI3_q = [];
Res.dynErr_PFI4_q = [];

Res.stErr_PFI1_H = [];
Res.stErr_PFI2_H = [];
Res.stErr_PFI3_H = [];
Res.stErr_PFI4_H = [];

Res.dynErr_PFI1_H = [];
Res.dynErr_PFI2_H = [];
Res.dynErr_PFI3_H = [];
Res.dynErr_PFI4_H = [];

%--------------------------------------------------------------------------
% Simulation
%--------------------------------------------------------------------------
for i=1:Res.nsamples
    % Euler Equ. Errors____________________________________________________
    for order_ct=1:M
        Poly = load(strcat(folder,'/Poly',num2str(order_ct-1),'_PFI.mat'));
        Sol = load(strcat(folder,'/Sol',num2str(order_ct-1),'_PFI.mat'));
        
        % Algorithm-implied policies
        pdf_aux = (squeeze(pdf(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)])';
        A(order_ct,i) = StaticParams.aGrid*pdf_aux';
        sd(order_ct,i) = sqrt(StaticParams.aGrid.^2*pdf_aux');
        skew(order_ct,i) = ((StaticParams.aGrid-squeeze(A(order_ct,i))).^3*pdf_aux')...
                            ./squeeze(sd(order_ct,i))^3;
        kurt(order_ct,i) = ((StaticParams.aGrid-squeeze(A(order_ct,i))).^4*pdf_aux')...
                            ./squeeze(sd(order_ct,i))^4;
        gini(order_ct,i) = sum(sum((pdf_aux'*pdf_aux).*abs(repmat(StaticParams.aGrid',[1,length(pdf_aux)])-repmat(StaticParams.aGrid,[length(pdf_aux),1])),1),2)...
                           /(2.*StaticParams.aGrid*pdf_aux');
        h_approx = interp1(StaticParams.aGrid_pol,squeeze(...
                  [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))')])',StaticParams.aGrid)';
        H(order_ct,i) = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]*sum(h_approx.*squeeze(pdf(order_ct,:,:)),2);
        b_approx = interp1(StaticParams.aGrid_pol,squeeze(...
                  [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.b(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.b(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))')])',StaticParams.aGrid)';
        b_aux = b_approx.*squeeze(pdf(order_ct,:,:));
        Borrow(order_ct,i) = StaticParams.ur(Res.z_ag(i)+1)*sum(b_aux(1,b_approx(1,:)<=0))+StaticParams.er(Res.z_ag(i)+1)*sum(b_aux(2,b_approx(2,:)<=0));
        Lend(order_ct,i) = StaticParams.ur(Res.z_ag(i)+1)*sum(b_aux(1,b_approx(1,:)>0))+StaticParams.er(Res.z_ag(i)+1)*sum(b_aux(2,b_approx(2,:)>0));

        a_prime = interp1(StaticParams.aGrid_pol,squeeze(...
                  [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))')])',Res.aGrid)';
        c_prime = interp1(StaticParams.aGrid_pol,squeeze(...
                  [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))')])',Res.aGrid)';   
        q(order_ct,i) = getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.q_pr(:,Res.z_ag(i)+1)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1))));
        r(order_ct,i) = getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,Res.z_ag(i)+1)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1))));
        if i==1
            ab_0 = A(order_ct,1);
            ag_0 = A(order_ct,1);
        else
            ab_0 = max(StaticParams.a_min,ab(order_ct,i-1));
            ag_0 = max(StaticParams.a_min,ag(order_ct,i-1));          
        end
        ab(order_ct,i) = interp1(StaticParams.aGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))'))',ab_0)';
        ag(order_ct,i) = interp1(StaticParams.aGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))'))',ag_0)';
        hb(order_ct,i) = interp1(StaticParams.aGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))'))',ab_0)';
        hg(order_ct,i) = interp1(StaticParams.aGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))'))',ag_0)';
        bb(order_ct,i) = ab(order_ct,i)-A(order_ct,i)*hb(order_ct,i);
        bg(order_ct,i) = ag(order_ct,i)-A(order_ct,i)*hg(order_ct,i);
        cb(order_ct,i) = interp1(StaticParams.aGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))'))',ab_0)';   
        cg(order_ct,i) = interp1(StaticParams.aGrid_pol,squeeze(...
                  getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))'))',ag_0)';   
        
        % standard Euler policies
        euler_weights = calcWeights(squeeze(cdf(order_ct,:,:)),Poly,StaticParams,Res.z_ag(i),...
                            interp1(Res.aGrid,sort(a_prime,2)',StaticParams.aGrid)');
        euler_expec = @(grid) Aggr(interp1(StaticParams.aGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')',...
                               interp1(StaticParams.aGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')')...
                      .^(-StaticParams.gamma)...
                      .*Aggr_derC(interp1(StaticParams.aGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')',...
                               interp1(StaticParams.aGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')');
        euler_r_pr = [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,1)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,1)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,2)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,2)),Sol.idx,squeeze(euler_weights(2,:))')];
        euler_c_prime = [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*((1+repmat(euler_r_pr,[1,size(a_prime,2)])).*euler_expec(a_prime(1,:))));...
                         (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*((1+repmat(euler_r_pr,[1,size(a_prime,2)])).*euler_expec(a_prime(2,:))))].^(-1/StaticParams.gamma);
        a_approx = interp1(StaticParams.aGrid_pol,squeeze(...
                  [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))');...
                   getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(weights(order_ct,1:(Poly.phi_N+1)))')])',StaticParams.aGrid)';
        euler_q = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]*sum(a_approx.*squeeze(pdf(order_ct,:,:)),2);
        
        % standard error
        Err = log10(abs(1-euler_c_prime./c_prime));%(c_prime-euler_c_prime)./max(eps,euler_c_prime);
        Err_algoGrid = log10(abs(1-euler_c_prime(:,algo_idx)./c_prime(:,algo_idx)));%(c_prime(:,algo_idx)-euler_c_prime(:,algo_idx))./max(eps,euler_c_prime(:,algo_idx));
        Err = Err(logical(squeeze(idx_posProb(order_ct,:,:))));
        Err_algoGrid = Err_algoGrid(logical(squeeze(idx_posProb(order_ct,:,algo_idx))));
        Err_q = log10(abs(1-euler_q./q(order_ct,i)));%(pr(order_ct,i)-euler_p)./max(eps,euler_p);
        Err_H = log10(abs(1-1./H(order_ct,i)));
        switch order_ct
            case 1
            Res.stErr_PFI1 = [Res.stErr_PFI1;Err];
            Res.stErr_PFI1_algoGrid = [Res.stErr_PFI1_algoGrid;Err_algoGrid];
            Res.stErr_PFI1_q = [Res.stErr_PFI1_q;Err_q];
            Res.stErr_PFI1_H = [Res.stErr_PFI1_H;Err_H];
            case 2
            Res.stErr_PFI2 = [Res.stErr_PFI2;Err];
            Res.stErr_PFI2_algoGrid = [Res.stErr_PFI2_algoGrid;Err_algoGrid];
            Res.stErr_PFI2_q = [Res.stErr_PFI2_q;Err_q];
            Res.stErr_PFI2_H = [Res.stErr_PFI2_H;Err_H];
            case 3
            Res.stErr_PFI3 = [Res.stErr_PFI3;Err];
            Res.stErr_PFI3_algoGrid = [Res.stErr_PFI3_algoGrid;Err_algoGrid];
            Res.stErr_PFI3_q = [Res.stErr_PFI3_q;Err_q];
            Res.stErr_PFI3_H = [Res.stErr_PFI3_H;Err_H];
            case 4
            Res.stErr_PFI4 = [Res.stErr_PFI4;Err];
            Res.stErr_PFI4_algoGrid = [Res.stErr_PFI4_algoGrid;Err_algoGrid];
            Res.stErr_PFI4_q = [Res.stErr_PFI4_q;Err_q];
            Res.stErr_PFI4_H = [Res.stErr_PFI4_H;Err_H];
        end
        % update distribution
        if i<Res.nsamples
            pdf(order_ct,:,:) = calcEoPDistr(interp1(Res.aGrid,sort(a_prime,2)',StaticParams.aGrid)',...
                           squeeze(cdf(order_ct,:,:)),StaticParams.aGrid);
            pdf(order_ct,:,:) = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                           .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                           *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(pdf,4)])...
                           .*squeeze(pdf(order_ct,:,:)));
            idx = max(find(squeeze(pdf(order_ct,1,:))>eps,1,'last'),find(squeeze(pdf(order_ct,2,:))>eps,1,'last'));
            idx_posProb(order_ct,:,(Res.aGrid<=StaticParams.aGrid(idx))) = 1;
            pdf = pdf.*(pdf>eps);
            pdf = pdf./repmat(sum(pdf,3),[1,1,size(pdf,3)]);
            cdf = cummax(min(1,cumsum(pdf,3)),4);
            cdf(order_ct,1,cdf(order_ct,1,:)==max(cdf(order_ct,1,:))) = 1;
            cdf(order_ct,2,cdf(order_ct,2,:)==max(cdf(order_ct,2,:))) = 1;
            weights(order_ct,1:(Poly.phi_N+1))...
            = calcWeights(squeeze(cdf(order_ct,:,:)),Poly,StaticParams,Res.z_ag(i+1));
        end
        
        % dynamic Euler policies and distributions
        dyn_euler_A(order_ct,i) = StaticParams.aGrid*squeeze(dyn_pdf(order_ct,:,:))'*[StaticParams.ur(Res.z_ag(i)+1);StaticParams.er(Res.z_ag(i)+1)];
        a_prime_next = [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,2*Res.z_ag(i)+1,:)),Sol.idx,squeeze(dyn_weights(order_ct,1:(Poly.phi_N+1)))');...
                        getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.a_prime(:,2*Res.z_ag(i)+2,:)),Sol.idx,squeeze(dyn_weights(order_ct,1:(Poly.phi_N+1)))')];  
        euler_weights = calcWeights(squeeze(dyn_cdf(order_ct,:,:)),Poly,StaticParams,Res.z_ag(i),...
                        interp1(StaticParams.aGrid_pol,sort(a_prime_next,2)',StaticParams.aGrid)');
        euler_expec = @(grid) Aggr(interp1(StaticParams.aGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')',...
                               interp1(StaticParams.aGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')')...
                      .^(-StaticParams.gamma)...
                      .*Aggr_derC(interp1(StaticParams.aGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')',...
                               interp1(StaticParams.aGrid_pol,...
                               [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,1,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,2,:)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,3,:)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.h(:,4,:)),Sol.idx,squeeze(euler_weights(2,:))')]',...
                               grid,'linear','extrap')');
        euler_r_pr = [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,1)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,1)),Sol.idx,squeeze(euler_weights(1,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,2)),Sol.idx,squeeze(euler_weights(2,:))');...
                                getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.r(:,2)),Sol.idx,squeeze(euler_weights(2,:))')];
        euler_q = 0.5*sum(euler_weights(:,1));
        dyn_euler_c_prime = [(StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+1,:)*((1+repmat(euler_r_pr,[1,size(a_prime,2)])).*euler_expec(a_prime(1,:))));...
                         (StaticParams.beta*StaticParams.P(2*Res.z_ag(i)+2,:)*((1+repmat(euler_r_pr,[1,size(a_prime,2)])).*euler_expec(a_prime(2,:))))].^(-1/StaticParams.gamma);
        dyn_euler_h_prime = F_h(r(order_ct,i),dyn_euler_A(order_ct,i),euler_q,dyn_euler_c_prime);
        dyn_euler_wealth = repmat(StaticParams.wage(2*Res.z_ag(i)+(1:2),1),[1,size(a_large,2)])+(1+r(order_ct,i)).*a_large;
        dyn_euler_a_prime = max(0,dyn_euler_wealth-dyn_euler_c_prime-F_uq(r(order_ct,i),dyn_euler_A(order_ct,i),euler_q).*dyn_euler_h_prime);
        a_prime = interp1(Res.aGrid,dyn_euler_a_prime',StaticParams.aGrid)';
        dyn_euler_q = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]*sum(a_prime.*squeeze(dyn_pdf(order_ct,:,:)),2);
        h_prime = interp1(Res.aGrid,dyn_euler_h_prime',StaticParams.aGrid)';
        dyn_euler_H(order_ct,i) = [StaticParams.ur(Res.z_ag(i)+1),StaticParams.er(Res.z_ag(i)+1)]*sum(h_prime.*squeeze(dyn_pdf(order_ct,:,:)),2);
        b_approx = interp1(Res.aGrid,(dyn_euler_a_prime-dyn_euler_q*dyn_euler_h_prime)',StaticParams.aGrid)';
        b_aux = b_approx.*squeeze(pdf(order_ct,:,:));
        dyn_euler_Borrow(order_ct,i) = StaticParams.ur(Res.z_ag(i)+1)*sum(b_aux(1,b_approx(1,:)<=0))+StaticParams.er(Res.z_ag(i)+1)*sum(b_aux(2,b_approx(2,:)<=0));
        dyn_euler_Lend(order_ct,i) = StaticParams.ur(Res.z_ag(i)+1)*sum(b_aux(1,b_approx(1,:)>0))+StaticParams.er(Res.z_ag(i)+1)*sum(b_aux(2,b_approx(2,:)>0));

        % dynamic error
        Err = log10(abs(1-dyn_euler_c_prime./c_prime));%(c_prime-dyn_euler_c_prime)./max(eps,dyn_euler_c_prime);
        Err_algoGrid = log10(abs(1-dyn_euler_c_prime(:,algo_idx)./c_prime(:,algo_idx)));%(c_prime(:,algo_idx)-dyn_euler_c_prime(:,algo_idx))./max(eps,dyn_euler_c_prime(:,algo_idx));
        Err = Err(logical(squeeze(dyn_idx_posProb(order_ct,:,:))));
        Err_algoGrid = Err_algoGrid(logical(squeeze(dyn_idx_posProb(order_ct,:,algo_idx))));
        Err_q = log10(abs(1-dyn_euler_q./q(order_ct,i)));%(pr(order_ct,i)-dyn_euler_p)./max(eps,dyn_euler_p);
        Err_H = log10(abs(1-1./dyn_euler_H(order_ct,i)));
        switch order_ct
            case 1
            Res.dynErr_PFI1 = [Res.dynErr_PFI1;Err];
            Res.dynErr_PFI1_algoGrid = [Res.dynErr_PFI1_algoGrid;Err_algoGrid];
            Res.dynErr_PFI1_q = [Res.dynErr_PFI1_q;Err_q];
            Res.dynErr_PFI1_H = [Res.dynErr_PFI1_H;Err_H];
            case 2
            Res.dynErr_PFI2 = [Res.dynErr_PFI2;Err];
            Res.dynErr_PFI2_algoGrid = [Res.dynErr_PFI2_algoGrid;Err_algoGrid];
            Res.dynErr_PFI2_q = [Res.dynErr_PFI2_q;Err_q];
            Res.dynErr_PFI2_H = [Res.dynErr_PFI2_H;Err_H];
            case 3
            Res.dynErr_PFI3 = [Res.dynErr_PFI3;Err];
            Res.dynErr_PFI3_algoGrid = [Res.dynErr_PFI3_algoGrid;Err_algoGrid];
            Res.dynErr_PFI3_q = [Res.dynErr_PFI3_q;Err_q];
            Res.dynErr_PFI3_H = [Res.dynErr_PFI3_H;Err_H];
            case 4
            Res.dynErr_PFI4 = [Res.dynErr_PFI4;Err];
            Res.dynErr_PFI4_algoGrid = [Res.dynErr_PFI4_algoGrid;Err_algoGrid];
            Res.dynErr_PFI4_q = [Res.dynErr_PFI4_q;Err_q];
            Res.dynErr_PFI4_H = [Res.dynErr_PFI4_H;Err_H];
        end
        % update distribution
        if i<Res.nsamples
            dyn_pdf(order_ct,:,:) = calcEoPDistr(interp1(Res.aGrid,sort(dyn_euler_a_prime,2)',StaticParams.aGrid)',...
                           squeeze(dyn_cdf(order_ct,:,:)),StaticParams.aGrid);
            dyn_pdf(order_ct,:,:) = (repmat(1./(StaticParams.p(2*Res.z_ag(i)+(1:2))'*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2)))',[1,2])...
                           .*StaticParams.P(2*Res.z_ag(i)+(1:2),2*Res.z_ag(i+1)+(1:2))')...
                           *(repmat(StaticParams.p(2*Res.z_ag(i)+(1:2)),[1,size(dyn_pdf,4)])...
                           .*squeeze(dyn_pdf(order_ct,:,:)));
            idx = max(find(squeeze(dyn_pdf(order_ct,1,:))>eps,1,'last'),find(squeeze(dyn_pdf(order_ct,2,:))>eps,1,'last'));
            dyn_idx_posProb(order_ct,:,(Res.aGrid<=StaticParams.aGrid(idx))) = 1;
            dyn_pdf = dyn_pdf.*(dyn_pdf>eps);
            dyn_pdf = dyn_pdf./repmat(sum(dyn_pdf,3),[1,1,size(dyn_pdf,3)]);
            dyn_cdf = cummax(min(1,cumsum(dyn_pdf,3)),3);
            dyn_cdf(order_ct,1,dyn_cdf(order_ct,1,:)==max(dyn_cdf(order_ct,1,:))) = 1;
            dyn_cdf(order_ct,2,dyn_cdf(order_ct,2,:)==max(dyn_cdf(order_ct,2,:))) = 1;
            dyn_weights(order_ct,1:(Poly.phi_N+1))...
            = calcWeights(squeeze(dyn_cdf(order_ct,:,:)),Poly,StaticParams,Res.z_ag(i+1));
        end
    end
    
    if i<Res.nsamples    
        if Res.z_ag(i)==0
            fill(x,y,'r','EdgeColor','none')
        else
            fill(x,y,'g','EdgeColor','none')
        end
        set(gca,'nextplot','add');
        pdf_aux = [squeeze(pdf(M,1,:))',squeeze(pdf(M,1,2:end))'];
        pdf_aux = pdf_aux(idx_pdf);
        plot(grid,pdf_aux,'m','LineWidth',2);
        set(gca,'nextplot','add');
        pdf_aux = [squeeze(pdf(M,2,:))',squeeze(pdf(M,2,2:end))'];
        pdf_aux = pdf_aux(idx_pdf);
        plot(grid,pdf_aux,'b','LineWidth',2);
        set(gca,'nextplot','replacechildren');
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
close(v)

Res.ab = ab;
Res.ag = ag;
Res.hb = hb;
Res.hg = hg;
Res.bb = bb;
Res.bg = bg;
Res.cb = cb;
Res.cg = cg;

Res.q = q;
Res.q_mean = mean(Res.q,2);
Res.q_mean0 = mean(Res.q(:,(Res.z_ag==0)),2);
Res.q_mean1 = mean(Res.q(:,(Res.z_ag==1)),2);
Res.q_sd = std(Res.q,[],2);
Res.q_sd0 = std(Res.q(:,(Res.z_ag==0)),[],2);
Res.q_sd1 = std(Res.q(:,(Res.z_ag==1)),[],2);
Res.q_skew = skewness(Res.q,[],2);
Res.q_skew0 = skewness(Res.q(:,(Res.z_ag==0)),[],2);
Res.q_skew1 = skewness(Res.q(:,(Res.z_ag==1)),[],2);
Res.q_kurt = kurtosis(Res.q,[],2);
Res.q_kurt0 = kurtosis(Res.q(:,(Res.z_ag==0)),[],2);
Res.q_kurt1 = kurtosis(Res.q(:,(Res.z_ag==1)),[],2);
Res.r = r;
Res.r_mean = mean(Res.r,2);
Res.r_mean0 = mean(Res.r(:,(Res.z_ag==0)),2);
Res.r_mean1 = mean(Res.r(:,(Res.z_ag==1)),2);
Res.r_sd = std(Res.r,[],2);
Res.r_sd0 = std(Res.r(:,(Res.z_ag==0)),[],2);
Res.r_sd1 = std(Res.r(:,(Res.z_ag==1)),[],2);
Res.r_skew = skewness(Res.r,[],2);
Res.r_skew0 = skewness(Res.r(:,(Res.z_ag==0)),[],2);
Res.r_skew1 = skewness(Res.r(:,(Res.z_ag==1)),[],2);
Res.r_kurt = kurtosis(Res.r,[],2);
Res.r_kurt0 = kurtosis(Res.r(:,(Res.z_ag==0)),[],2);
Res.r_kurt1 = kurtosis(Res.r(:,(Res.z_ag==1)),[],2);
Res.A = A;
Res.sd = sd;
Res.skew = skew;
Res.kurt = kurt;
Res.gini = gini;
Res.Borrow = Borrow;
Res.Lend = Lend;
Res.dyn_euler_A = dyn_euler_A;

%--------------------------------------------------------------------------
% Compute the Euler equation error statistics
%--------------------------------------------------------------------------
Res.stErr_PFI1 = Res.stErr_PFI1(~isnan(Res.stErr_PFI1));
Res.stErr_PFI2 = Res.stErr_PFI2(~isnan(Res.stErr_PFI2));
Res.stErr_PFI3 = Res.stErr_PFI3(~isnan(Res.stErr_PFI3));
Res.stErr_PFI4 = Res.stErr_PFI4(~isnan(Res.stErr_PFI4));
Res.dynErr_PFI1 = Res.dynErr_PFI1(~isnan(Res.dynErr_PFI1));
Res.dynErr_PFI2 = Res.dynErr_PFI2(~isnan(Res.dynErr_PFI2));
Res.dynErr_PFI3 = Res.dynErr_PFI3(~isnan(Res.dynErr_PFI3));
Res.dynErr_PFI4 = Res.dynErr_PFI4(~isnan(Res.dynErr_PFI4));
Res.stErr_PFI1_algoGrid = Res.stErr_PFI1_algoGrid(~isnan(Res.stErr_PFI1_algoGrid));
Res.stErr_PFI2_algoGrid = Res.stErr_PFI2_algoGrid(~isnan(Res.stErr_PFI2_algoGrid));
Res.stErr_PFI3_algoGrid = Res.stErr_PFI3_algoGrid(~isnan(Res.stErr_PFI3_algoGrid));
Res.stErr_PFI4_algoGrid = Res.stErr_PFI4_algoGrid(~isnan(Res.stErr_PFI4_algoGrid));
Res.dynErr_PFI1_algoGrid = Res.dynErr_PFI1_algoGrid(~isnan(Res.dynErr_PFI1_algoGrid));
Res.dynErr_PFI2_algoGrid = Res.dynErr_PFI2_algoGrid(~isnan(Res.dynErr_PFI2_algoGrid));
Res.dynErr_PFI3_algoGrid = Res.dynErr_PFI3_algoGrid(~isnan(Res.dynErr_PFI3_algoGrid));
Res.dynErr_PFI4_algoGrid = Res.dynErr_PFI4_algoGrid(~isnan(Res.dynErr_PFI4_algoGrid));
Res.median_stErr_PFI = [median((Res.stErr_PFI1));median((Res.stErr_PFI2));median((Res.stErr_PFI3));median((Res.stErr_PFI4))];
Res.mean_stErr_PFI = [mean((Res.stErr_PFI1));mean((Res.stErr_PFI2));mean((Res.stErr_PFI3));mean((Res.stErr_PFI4))];
Res.min_stErr_PFI = [min((Res.stErr_PFI1));min((Res.stErr_PFI2));min((Res.stErr_PFI3));min((Res.stErr_PFI4))];
Res.max_stErr_PFI = [max((Res.stErr_PFI1));max((Res.stErr_PFI2));max((Res.stErr_PFI3));max((Res.stErr_PFI4))];
Res.median_dynErr_PFI = [median((Res.dynErr_PFI1));median((Res.dynErr_PFI2));median((Res.dynErr_PFI3));median((Res.dynErr_PFI4))];
Res.mean_dynErr_PFI = [mean((Res.dynErr_PFI1));mean((Res.dynErr_PFI2));mean((Res.dynErr_PFI3));mean((Res.dynErr_PFI4))];
Res.min_dynErr_PFI = [min((Res.dynErr_PFI1));min((Res.dynErr_PFI2));min((Res.dynErr_PFI3));min((Res.dynErr_PFI4))];
Res.max_dynErr_PFI = [max((Res.dynErr_PFI1));max((Res.dynErr_PFI2));max((Res.dynErr_PFI3));max((Res.dynErr_PFI4))];

Res.median_stErr_PFI_algoGrid = [median((Res.stErr_PFI1_algoGrid));median((Res.stErr_PFI2_algoGrid));median((Res.stErr_PFI3_algoGrid));median((Res.stErr_PFI4_algoGrid))];
Res.mean_stErr_PFI_algoGrid = [mean((Res.stErr_PFI1_algoGrid));mean((Res.stErr_PFI2_algoGrid));mean((Res.stErr_PFI3_algoGrid));mean((Res.stErr_PFI4_algoGrid))];
Res.min_stErr_PFI_algoGrid = [min((Res.stErr_PFI1_algoGrid));min((Res.stErr_PFI2_algoGrid));min((Res.stErr_PFI3_algoGrid));min((Res.stErr_PFI4_algoGrid))];
Res.max_stErr_PFI_algoGrid = [max((Res.stErr_PFI1_algoGrid));max((Res.stErr_PFI2_algoGrid));max((Res.stErr_PFI3_algoGrid));max((Res.stErr_PFI4_algoGrid))];
Res.median_dynErr_PFI_algoGrid = [median((Res.dynErr_PFI1_algoGrid));median((Res.dynErr_PFI2_algoGrid));median((Res.dynErr_PFI3_algoGrid));median((Res.dynErr_PFI4_algoGrid))];
Res.mean_dynErr_PFI_algoGrid = [mean((Res.dynErr_PFI1_algoGrid));mean((Res.dynErr_PFI2_algoGrid));mean((Res.dynErr_PFI3_algoGrid));mean((Res.dynErr_PFI4_algoGrid))];
Res.min_dynErr_PFI_algoGrid = [min((Res.dynErr_PFI1_algoGrid));min((Res.dynErr_PFI2_algoGrid));min((Res.dynErr_PFI3_algoGrid));min((Res.dynErr_PFI4_algoGrid))];
Res.max_dynErr_PFI_algoGrid = [max((Res.dynErr_PFI1_algoGrid));max((Res.dynErr_PFI2_algoGrid));max((Res.dynErr_PFI3_algoGrid));max((Res.dynErr_PFI4_algoGrid))];

Res.median_stErr_PFI_q = [median((Res.stErr_PFI1_q));median((Res.stErr_PFI2_q));median((Res.stErr_PFI3_q));median((Res.stErr_PFI4_q))];
Res.mean_stErr_PFI_q = [mean((Res.stErr_PFI1_q));mean((Res.stErr_PFI2_q));mean((Res.stErr_PFI3_q));mean((Res.stErr_PFI4_q))];
Res.min_stErr_PFI_q = [min((Res.stErr_PFI1_q));min((Res.stErr_PFI2_q));min((Res.stErr_PFI3_q));min((Res.stErr_PFI4_q))];
Res.max_stErr_PFI_q = [max((Res.stErr_PFI1_q));max((Res.stErr_PFI2_q));max((Res.stErr_PFI3_q));max((Res.stErr_PFI4_q))];
Res.median_dynErr_PFI_q = [median((Res.dynErr_PFI1_q));median((Res.dynErr_PFI2_q));median((Res.dynErr_PFI3_q));median((Res.dynErr_PFI4_q))];
Res.mean_dynErr_PFI_q = [mean((Res.dynErr_PFI1_q));mean((Res.dynErr_PFI2_q));mean((Res.dynErr_PFI3_q));mean((Res.dynErr_PFI4_q))];
Res.min_dynErr_PFI_q = [min((Res.dynErr_PFI1_q));min((Res.dynErr_PFI2_q));min((Res.dynErr_PFI3_q));min((Res.dynErr_PFI4_q))];
Res.max_dynErr_PFI_q = [max((Res.dynErr_PFI1_q));max((Res.dynErr_PFI2_q));max((Res.dynErr_PFI3_q));max((Res.dynErr_PFI4_q))];

Res.median_stErr_PFI_H = [median((Res.stErr_PFI1_H));median((Res.stErr_PFI2_H));median((Res.stErr_PFI3_H));median((Res.stErr_PFI4_H))];
Res.mean_stErr_PFI_H = [mean((Res.stErr_PFI1_H));mean((Res.stErr_PFI2_H));mean((Res.stErr_PFI3_H));mean((Res.stErr_PFI4_H))];
Res.min_stErr_PFI_H = [min((Res.stErr_PFI1_H));min((Res.stErr_PFI2_H));min((Res.stErr_PFI3_H));min((Res.stErr_PFI4_H))];
Res.max_stErr_PFI_H = [max((Res.stErr_PFI1_H));max((Res.stErr_PFI2_H));max((Res.stErr_PFI3_H));max((Res.stErr_PFI4_H))];
Res.median_dynErr_PFI_H = [median((Res.dynErr_PFI1_H));median((Res.dynErr_PFI2_H));median((Res.dynErr_PFI3_H));median((Res.dynErr_PFI4_H))];
Res.mean_dynErr_PFI_H = [mean((Res.dynErr_PFI1_H));mean((Res.dynErr_PFI2_H));mean((Res.dynErr_PFI3_H));mean((Res.dynErr_PFI4_H))];
Res.min_dynErr_PFI_H = [min((Res.dynErr_PFI1_H));min((Res.dynErr_PFI2_H));min((Res.dynErr_PFI3_H));min((Res.dynErr_PFI4_H))];
Res.max_dynErr_PFI_H = [max((Res.dynErr_PFI1_H));max((Res.dynErr_PFI2_H));max((Res.dynErr_PFI3_H));max((Res.dynErr_PFI4_H))];

save(strcat(folder,'/Res.mat'),'-struct','Res');
end