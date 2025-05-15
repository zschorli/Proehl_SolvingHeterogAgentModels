% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file produces the plots and tables in the paper.
%__________________________________________________________________________
clc;
clear all;

%--------------------------------------------------------------------------
% Boxplot of law of motion errors for the baseline
%--------------------------------------------------------------------------
case_nr=1;
folder = strcat('res_case',num2str(case_nr));
Res = load(strcat(folder,'/Res.mat'));
figure(2)
subplot(1,2,1)
a1=[Res.LoMErr_PFI0';Res.LoMErr_PFI1';Res.LoMErr_PFI2';Res.LoMErr_PFI3'];%Res.LoMErr_PFI4];
b1=[zeros(length(Res.LoMErr_PFI0),1);ones(length(Res.LoMErr_PFI1),1);...
   2*ones(length(Res.LoMErr_PFI2),1);3*ones(length(Res.LoMErr_PFI3),1)];%13*ones(length(Res.LoMErr_PFI4),1)];
boxplot(a1,b1,...
         'Labels',{'PFI0','PFI1','PFI2','PFI3'},...
         'OutlierSize',1,'Symbol','r.');%'Whisker',4
axis([0.5 4.5 -7.5 -2])
title('Panel A')

subplot(1,2,2)
a1=[Res.LoMErrL_PFI0';Res.LoMErrL_PFI1';Res.LoMErrL_PFI2';Res.LoMErrL_PFI3'];%Res.LoMErr_PFI4];
b1=[zeros(length(Res.LoMErrL_PFI0),1);ones(length(Res.LoMErrL_PFI1),1);...
   2*ones(length(Res.LoMErrL_PFI2),1);3*ones(length(Res.LoMErrL_PFI3),1)];%13*ones(length(Res.LoMErr_PFI4),1)];
boxplot(a1,b1,...
         'Labels',{'PFI0','PFI1','PFI2','PFI3'},...
         'OutlierSize',1,'Symbol','r.');%'Whisker',4
axis([0.5 4.5 -7.5 -2])
title('Panel B')

StaticParams = load(strcat(folder,'/StaticParams.mat'));
Sol0 = load(strcat(folder,'/Sol0_PFI.mat'));
Sol1 = load(strcat(folder,'/Sol1_PFI.mat'));
Sol2 = load(strcat(folder,'/Sol2_PFI.mat'));
Sol3 = load(strcat(folder,'/Sol3_PFI.mat'));

p=[sum(StaticParams.p([1,3])),sum(StaticParams.p([2,4]))];

figure(3)
hold on
plot(StaticParams.kGrid,p*Sol0.pdf,'r','LineWidth',1);
plot(StaticParams.kGrid,p*Sol1.pdf,'c','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol3.pdf,'m--','LineWidth',1.5);
hold off 
axis([0 25 0 0.004])
legend('PFI 0^{th} order','PFI 1^{st} order','PFI 2^{nd} order','PFI 3^{rd} order')
xlabel('individual capital k');
ylabel('probability');

%--------------------------------------------------------------------------
% Full Comparison of Stationary Distributions
%--------------------------------------------------------------------------
file_loc = 'C:\Users\Elisabeth Proehl\Documents\GitHub\Proehl_SolvingHeterogAgentModels\Code\Proehl_GrowthModel\';
StaticParams1 = load(strcat(file_loc,'res_case',num2str(1),'/StaticParams.mat'));
StaticParams2 = load(strcat(file_loc,'res_case',num2str(2),'/StaticParams.mat'));
StaticParams3 = load(strcat(file_loc,'res_case',num2str(3),'/StaticParams.mat'));
Sol1 = load(strcat(file_loc,'res_case',num2str(1),'/Sol3_PPA.mat'));
Sol2 = load(strcat(file_loc,'res_case',num2str(2),'/Sol3_PPA.mat'));
Sol3 = load(strcat(file_loc,'res_case',num2str(3),'/Sol3_PPA.mat'));
compSol = load(strcat(file_loc,'compareDistr_NoAggShock/compareSol.mat'));
p=[sum(StaticParams1.p([1,3])),sum(StaticParams1.p([2,4]))];

% endog. labor, progressive taxation
Sol1eL = load(strcat('res_case1','/Sol3_PFI.mat'));
Sol3eL = load(strcat('res_case3','/Sol3_PFI.mat'));
StaticParams1eL = load(strcat('res_case1','/StaticParams.mat'));
p1=[sum(StaticParams1eL.p([1,3])),sum(StaticParams1eL.p([2,4]))];

% endog. labor, linear taxation
Sol2eL = load(strcat('res_case2','/Sol3_PFI.mat'));
Sol4eL = load(strcat('res_case4','/Sol3_PFI.mat'));
StaticParams2eL = load(strcat('res_case2','/StaticParams.mat'));
p2=[sum(StaticParams2eL.p([1,3])),sum(StaticParams2eL.p([2,4]))];

figure(4)
subplot(1,3,1)
pdf = compSol.pdf_b;
pdf(:,2:end) = pdf(:,2:end)./repmat(StaticParams1.kGrid(2:end)-StaticParams1.kGrid(1:end-1),[size(pdf,1),1]); 
hold on
plot(StaticParams1.kGrid,pdf(1,:));
plot(StaticParams1.kGrid,pdf(4,:),'b','LineWidth',1.5)
plot(StaticParams1.kGrid,pdf(9,:),'r','LineWidth',1.5);
plot(StaticParams1.kGrid,pdf(14,:),'m','LineWidth',1.5);
plot(StaticParams1.kGrid,pdf(19,:));
plot(compSol.fullInsurance_K_b*ones(1,2),[0,1]','c');
hold off
axis([0 120 0 0.03])
names = cell(1,length(compSol.muGrid));
for i=[1,4,9,14,19]
    names{i} = strcat('\nu=',num2str(compSol.muGrid(i)*100),'%');
end
legend(names{[1,4,9,14,19]},'full UI')
xlabel('individual capital k');
ylabel('probability');
title('Panel A')

subplot(1,3,2)
hold on 
pdf = p*Sol1.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams1.kGrid(2:end)-StaticParams1.kGrid(1:end-1)); 
plot(StaticParams1.kGrid,pdf,'b','LineWidth',1.5);
pdf = p*Sol2.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams2.kGrid(2:end)-StaticParams2.kGrid(1:end-1)); 
plot(StaticParams2.kGrid,pdf,'r','LineWidth',1.5);
pdf = p*Sol3.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams3.kGrid(2:end)-StaticParams3.kGrid(1:end-1)); 
plot(StaticParams3.kGrid,pdf,'m','LineWidth',1.5);
hold off
axis([0 120 0 0.03])
legend('\nu=15%','\nu=40%','\nu=65%')
xlabel('individual capital k');
ylabel('probability');
title('Panel B')

subplot(1,3,3)
hold on 
pdf = p2*Sol2eL.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams2eL.kGrid(2:end)-StaticParams2eL.kGrid(1:end-1)); 
plot(StaticParams2eL.kGrid,pdf,'b','LineWidth',1.5);
pdf = p2*Sol4eL.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams2eL.kGrid(2:end)-StaticParams2eL.kGrid(1:end-1)); 
plot(StaticParams2eL.kGrid,pdf,'r','LineWidth',1.5);
pdf = p1*Sol1eL.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams1eL.kGrid(2:end)-StaticParams1eL.kGrid(1:end-1)); 
plot(StaticParams1eL.kGrid,pdf,'b--','LineWidth',1);
pdf = p1*Sol3eL.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams1eL.kGrid(2:end)-StaticParams1eL.kGrid(1:end-1)); 
plot(StaticParams1eL.kGrid,pdf,'r--','LineWidth',1);
hold off
axis([0 30 0 0.11])
legend('\nu=15%, lin.','\nu=40%, lin.','\nu=15%, prog.','\nu=40%, prog.')
xlabel('individual capital k');
ylabel('probability');
title('Panel C')

% welfare
figure(5)
welf_b = log(max(eps,compSol.cGrid))*compSol.pdf_c_b';
welf_der_b = (welf_b(2:end)-welf_b(1:end-1))./(compSol.muGrid(2:end)-compSol.muGrid(1:end-1));
welf_der_b(2:end) = 0.5.*(welf_der_b(2:end)+welf_der_b(1:end-1));
welInterp_b = @(nu) interp1(compSol.muGrid,welf_b,nu);
FOCinterp_b = @(nu) interp1(compSol.muGrid(1:end-1),welf_der_b,nu);
optPol_b = fzero(FOCinterp_b,0.1);
welf_g = log(max(eps,compSol.cGrid))*compSol.pdf_c_g';
welf_der_g = (welf_g(2:end)-welf_g(1:end-1))./(compSol.muGrid(2:end)-compSol.muGrid(1:end-1));
welf_der_g(2:end) = 0.5.*(welf_der_g(2:end)+welf_der_g(1:end-1));
welInterp_g = @(nu) interp1(compSol.muGrid,welf_g,nu);
FOCinterp_g = @(nu) interp1(compSol.muGrid(1:end-1),welf_der_g,nu);
if max(welf_der_g)>0 && min(welf_der_g)<0 
    optPol_g = fzero(FOCinterp_g,0.1);
else
    optPol_g = 0;
end
subplot(1,2,1)
hold on
plot(compSol.muGrid(1:19),welf_b(1:19),'b')
h_b=plot(optPol_b,welInterp_b(optPol_b),'b*');
scatter(1,log(compSol.fullInsurance_C_b),[],'blue','filled')
plot(compSol.muGrid(1:19),welf_g(1:19),'r')
h_g=plot(optPol_g,welInterp_g(optPol_g),'r*');
scatter(0.96/0.9,log(compSol.fullInsurance_C_g),[],'red','filled')
hold off
legend([h_b,h_g],{strcat('\nu=',num2str(100*round(optPol_b,4)),'%'),strcat('\nu=',num2str(100*round(optPol_g,4)),'%')},'location','east')
xlabel('\nu');
title('Panel A')

subplot(1,2,2)
hold on
plot(compSol.muGrid(1:19),welf_der_b(1:19),'b')
plot(optPol_b,0,'b*')
plot(compSol.muGrid(1:19),welf_der_g(1:19),'r')
hold off
xlabel('\nu');
title('Panel B')

%--------------------------------------------------------------------------
% Policy Experiment
%--------------------------------------------------------------------------
for case_nr = 15:20
    folder = strcat('res_case',num2str(case_nr));
    StaticParams = load(strcat(folder,'/StaticParams.mat'));
    Sol = load(strcat(folder,'/Sol3_PFI.mat'));
    Poly = load(strcat(folder,'/Poly3_PFI.mat'));
    weights0 = calcWeights(Sol.cdf_cond(1:2,:),Poly,StaticParams,0);
    weights1 = calcWeights(Sol.cdf_cond(3:4,:),Poly,StaticParams,1);
    c_prime0 = max(eps,interp1(StaticParams.kGrid_pol,squeeze(...
                      [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,weights0);...
                       getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,weights0)])',StaticParams.kGrid)');
    c_prime1 = max(eps,interp1(StaticParams.kGrid_pol,squeeze(...
                      [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,weights1);...
                       getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,weights1)])',StaticParams.kGrid)');
    l_prime0 = min(1-eps,interp1(StaticParams.kGrid_pol,squeeze(...
                      [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,1,:)),Sol.idx,weights0);...
                       getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,2,:)),Sol.idx,weights0)])',StaticParams.kGrid)');
    l_prime1 = min(1-eps,interp1(StaticParams.kGrid_pol,squeeze(...
                      [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,3,:)),Sol.idx,weights1);...
                       getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.l(:,4,:)),Sol.idx,weights1)])',StaticParams.kGrid)');
    aggregator0 = c_prime0.^StaticParams.theta.*(1-l_prime0).^(1-StaticParams.theta);
    aggregator1 = c_prime1.^StaticParams.theta.*(1-l_prime1).^(1-StaticParams.theta);
    if StaticParams.gamma==1
        ut0 = log(aggregator0);
        ut1 = log(aggregator1);
    else
        ut0 = (aggregator0).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
        ut1 = (aggregator1).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
    end
    Sol.agC = [sum(c_prime0.*Sol.pdf_cond(1:2,:),2)',...
               sum(c_prime1.*Sol.pdf_cond(3:4,:),2)'];
    Sol.agWelfare = [sum(ut0.*Sol.pdf_cond(1:2,:),2)',...
                     sum(ut1.*Sol.pdf_cond(3:4,:),2)'];
    save(strcat(folder,'/Sol3_PFI.mat'),'-struct','Sol'); 
end 

Sol3 = load(strcat('res_case3','/Sol3_PFI.mat'));
Sol5 = load(strcat('res_case5','/Sol3_PFI.mat'));
Sol6 = load(strcat('res_case6','/Sol3_PFI.mat'));
Sol7 = load(strcat('res_case7','/Sol3_PFI.mat'));
% Sol8 = load(strcat('res_case8','/Sol3_PFI.mat'));
% Sol9 = load(strcat('res_case9','/Sol3_PFI.mat'));
% Sol10 = load(strcat('res_case10','/Sol3_PFI.mat'));
Sol11 = load(strcat('res_case11','/Sol3_PFI.mat'));
Sol12 = load(strcat('res_case12','/Sol3_PFI.mat'));

Sol15 = load(strcat('res_case15','/Sol3_PFI.mat'));
Sol16 = load(strcat('res_case16','/Sol3_PFI.mat'));
Sol17 = load(strcat('res_case17','/Sol3_PFI.mat'));
Sol18 = load(strcat('res_case18','/Sol3_PFI.mat'));
Sol19 = load(strcat('res_case19','/Sol3_PFI.mat'));
Sol20 = load(strcat('res_case20','/Sol3_PFI.mat'));

StaticParams = load(strcat('res_case3','/StaticParams.mat'));

pL = StaticParams.p;
pL([1,3])=0;
pL = pL./sum(pL);

N=6*4;
Output = zeros(15,N);

% aggregate capital
Output(1,1:4:N) = [(StaticParams.kGrid*(StaticParams.p'*Sol3.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol5.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol6.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol7.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol11.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol12.pdf_cond)')];
Output(2,1:2:N) = [(StaticParams.kGrid*(StaticParams.p(1:2)'*Sol3.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol3.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol5.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol5.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol6.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol6.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol7.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol7.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol11.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol11.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol12.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol12.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4))];
Output(3,:) = [(StaticParams.kGrid*Sol3.pdf_cond'),...
     (StaticParams.kGrid*Sol5.pdf_cond'),...
     (StaticParams.kGrid*Sol6.pdf_cond'),...
     (StaticParams.kGrid*Sol7.pdf_cond'),...
     (StaticParams.kGrid*Sol11.pdf_cond'),...
     (StaticParams.kGrid*Sol12.pdf_cond')];
     
% labor supply
Output(4,1:4:N) = [(StaticParams.lGrid*(pL'*Sol3.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol5.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol6.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol7.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol11.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol12.pdf_cond_l)')];
Output(5,1:2:N) = [(StaticParams.lGrid*(pL(1:2)'*Sol3.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol3.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol5.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol5.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol6.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol6.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol7.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol7.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol11.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol11.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol12.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol12.pdf_cond_l(3:4,:))')./sum(pL(3:4))];  

% output
Output(7,1:4:N) = StaticParams.output(0,Output(2,1:4:N),Output(5,1:4:N));
Output(7,3:4:N) = StaticParams.output(1,Output(3,2:4:N),Output(5,3:4:N));
aux = Output(7,1:2:N);
Output(6,1:4:N) = 0.5.*(aux(1:2:N/2)+aux(2:2:N/2));

% productivity: output per hour worked 
Output(9,:) = Output(7,:)./max(eps,Output(5,:));   
aux = Output(9,1:2:N);
Output(8,1:4:N) = 0.5.*(aux(1:2:N/2)+aux(2:2:N/2));

% consumption
Output(10,1:4:N) = [Sol3.agC*StaticParams.p,Sol5.agC*StaticParams.p,Sol6.agC*StaticParams.p,...
                    Sol7.agC*StaticParams.p,Sol11.agC*StaticParams.p,Sol12.agC*StaticParams.p];
Output(11,1:2:N) = 2.*[Sol3.agC(1:2)*StaticParams.p(1:2),Sol3.agC(3:4)*StaticParams.p(3:4),...
                     Sol5.agC(1:2)*StaticParams.p(1:2),Sol5.agC(3:4)*StaticParams.p(3:4),...
                     Sol6.agC(1:2)*StaticParams.p(1:2),Sol6.agC(3:4)*StaticParams.p(3:4),...
                     Sol7.agC(1:2)*StaticParams.p(1:2),Sol7.agC(3:4)*StaticParams.p(3:4),...
                     Sol11.agC(1:2)*StaticParams.p(1:2),Sol11.agC(3:4)*StaticParams.p(3:4),...
                     Sol12.agC(1:2)*StaticParams.p(1:2),Sol12.agC(3:4)*StaticParams.p(3:4)];
Output(12,:) = [Sol3.agC,Sol5.agC,Sol6.agC,...
                Sol7.agC,Sol11.agC,Sol12.agC];

% welfare
Output(13,1:4:N) = [Sol3.agWelfare*StaticParams.p,...
                    Sol5.agWelfare*StaticParams.p,...
                    Sol6.agWelfare*StaticParams.p,...
                    Sol7.agWelfare*StaticParams.p,...
                    Sol11.agWelfare*StaticParams.p,...
                    Sol12.agWelfare*StaticParams.p];
Output(14,1:2:N) = 2.*[Sol3.agWelfare(1:2)*StaticParams.p(1:2),Sol3.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol5.agWelfare(1:2)*StaticParams.p(1:2),Sol5.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol6.agWelfare(1:2)*StaticParams.p(1:2),Sol6.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol7.agWelfare(1:2)*StaticParams.p(1:2),Sol7.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol11.agWelfare(1:2)*StaticParams.p(1:2),Sol11.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol12.agWelfare(1:2)*StaticParams.p(1:2),Sol12.agWelfare(3:4)*StaticParams.p(3:4)];
Output(15,:) = [Sol3.agWelfare,Sol5.agWelfare,Sol6.agWelfare,...
                Sol7.agWelfare,Sol11.agWelfare,Sol12.agWelfare];
% derivatives
Output([1,4,6,8,10,13],6)=(Output([1,4,6,8,10,13],5)-Output([1,4,6,8,10,13],1))/norm(0.01*ones(1,2));
Output([1,4,6,8,10,13],10)=(Output([1,4,6,8,10,13],9)-Output([1,4,6,8,10,13],1))/norm(0.01*ones(1,2));
Output([1,4,6,8,10,13],18)=(Output([1,4,6,8,10,13],17)-Output([1,4,6,8,10,13],13))/norm(0.01*ones(1,2));
Output([1,4,6,8,10,13],22)=(Output([1,4,6,8,10,13],21)-Output([1,4,6,8,10,13],13))/norm(0.01*ones(1,2));

% Fit welfare function to 10 computed values, x grid: benfit in boom, 
% y grid: benefit in recession
x = zeros(10,1);
y = zeros(10,1);
w = zeros(10,1);
count = 0;
for case_nr=[1,3,5:12]
    count = count+1;
    folder = strcat('res_case',num2str(case_nr));
    StaticParams = load(strcat(folder,'/StaticParams.mat'));
    Sol = load(strcat(folder,'/Sol3_PFI.mat'));
    x(count) = StaticParams.mu;
    y(count) = StaticParams.mu+StaticParams.mu_extra;
    w(count) = Sol.agWelfare*StaticParams.p;
end
A = [ones(10,1),x,y,x.^2,y.^2,x.*y,x.^3,y.^3,x.*y.^2,x.^2.*y];
coeffs = linsolve(A,w);
welfareFct = @(x,y) coeffs(1)+coeffs(2).*x+coeffs(3).*y+coeffs(4).*x.^2 ...
                   +coeffs(5).*y.^2+coeffs(6).*x.*y+coeffs(7).*x.^3 ...
                   +coeffs(8).*y.^3+coeffs(9).*x.*y.^2+coeffs(10).*x.^2.*y; 
grid = linspace(0,0.5,51)';
[X,Y] = ndgrid(grid,grid);
Z = welfareFct(X,Y);
surf(X,Y,Z)