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
ct = num2cell(1:14);
for case_nr = [1,3,5:24,35:44]
    folder = strcat('res_case',num2str(case_nr));
    StaticParams = load(strcat(folder,'/StaticParams.mat'));
    switch case_nr
        case {ct{:},15,17,18,35,37,38}
            count = 3;
        otherwise
            count = 2;
    end
    Sol = load(strcat(folder,'/Sol',num2str(count),'_PFI.mat'));
    Poly = load(strcat(folder,'/Poly',num2str(count),'_PFI.mat'));
    weights0 = calcWeights(Sol.cdf_cond(1:2,:),Poly,StaticParams,0);
    weights1 = calcWeights(Sol.cdf_cond(3:4,:),Poly,StaticParams,1);
    c_prime0 = max(0,interp1(StaticParams.kGrid_pol,squeeze(...
                      [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,weights0);...
                       getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,weights0)])',StaticParams.kGrid)');
    c_prime0(c_prime0<eps^2) = 0;
    c_prime0(Sol.pdf_cond(1:2,:)==0) = max(eps,c_prime0(Sol.pdf_cond(1:2,:)==0));
    c_prime1 = max(0,interp1(StaticParams.kGrid_pol,squeeze(...
                      [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,weights1);...
                       getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,weights1)])',StaticParams.kGrid)');
    c_prime1(Sol.pdf_cond(3:4,:)==0) = max(eps,c_prime1(Sol.pdf_cond(3:4,:)==0));
    c_prime1(c_prime1<eps^2) = 0;
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
        ut0_b = log(max(eps^20,aggregator0));
        ut1_b = log(max(eps^20,aggregator1));
    else
        ut0 = (aggregator0).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
        ut1 = (aggregator1).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
        ut0_b = (max(eps^20,aggregator0)).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
        ut1_b = (max(eps^20,aggregator1)).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
    end
    Sol.agC = [sum(c_prime0.*Sol.pdf_cond(1:2,:),2)',...
               sum(c_prime1.*Sol.pdf_cond(3:4,:),2)'];
    Sol.agWelfare = [sum(ut0.*Sol.pdf_cond(1:2,:),2)',...
                     sum(ut1.*Sol.pdf_cond(3:4,:),2)'];
    Sol.agWelfare_b = [sum(ut0_b.*Sol.pdf_cond(1:2,:),2)',...
                       sum(ut1_b.*Sol.pdf_cond(3:4,:),2)'];
    save(strcat(folder,'/Sol',num2str(count),'_PFI.mat'),'-struct','Sol'); 
end 

Sol1 = load(strcat('res_case1','/Sol3_PFI.mat'));
Sol5 = load(strcat('res_case5','/Sol3_PFI.mat'));
Sol6 = load(strcat('res_case6','/Sol3_PFI.mat'));
Sol7 = load(strcat('res_case7','/Sol3_PFI.mat'));
Sol11 = load(strcat('res_case11','/Sol3_PFI.mat'));
Sol12 = load(strcat('res_case12','/Sol3_PFI.mat'));
Sol15 = load(strcat('res_case35','/Sol3_PFI.mat'));
Sol17 = load(strcat('res_case37','/Sol3_PFI.mat'));
Sol18 = load(strcat('res_case38','/Sol3_PFI.mat'));

StaticParams = load(strcat('res_case1','/StaticParams.mat'));
StaticParamsCC = load(strcat('res_case35','/StaticParams.mat'));

pL = StaticParams.p;
pL([1,3])=0;
pL = pL./sum(pL);

N=9*4;
Output = zeros(15,N);

% aggregate capital
Output(1,1:4:N) = [(StaticParams.kGrid*(StaticParams.p'*Sol1.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol5.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol6.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol7.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol11.pdf_cond)'),...
     (StaticParams.kGrid*(StaticParams.p'*Sol12.pdf_cond)'),...
     (StaticParamsCC.kGrid*(StaticParamsCC.p'*Sol15.pdf_cond)'),...
     (StaticParamsCC.kGrid*(StaticParamsCC.p'*Sol17.pdf_cond)'),...
     (StaticParamsCC.kGrid*(StaticParamsCC.p'*Sol18.pdf_cond)')];
Output(2,1:2:N) = [(StaticParams.kGrid*(StaticParams.p(1:2)'*Sol1.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol1.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol5.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol5.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol6.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol6.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol7.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol7.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol11.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol11.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParams.kGrid*(StaticParams.p(1:2)'*Sol12.pdf_cond(1:2,:))')./sum(StaticParams.p(1:2)),...
     (StaticParams.kGrid*(StaticParams.p(3:4)'*Sol12.pdf_cond(3:4,:))')./sum(StaticParams.p(3:4)),...
     (StaticParamsCC.kGrid*(StaticParamsCC.p(1:2)'*Sol15.pdf_cond(1:2,:))')./sum(StaticParamsCC.p(1:2)),...
     (StaticParamsCC.kGrid*(StaticParamsCC.p(3:4)'*Sol15.pdf_cond(3:4,:))')./sum(StaticParamsCC.p(3:4)),...
     (StaticParamsCC.kGrid*(StaticParamsCC.p(1:2)'*Sol17.pdf_cond(1:2,:))')./sum(StaticParamsCC.p(1:2)),...
     (StaticParamsCC.kGrid*(StaticParamsCC.p(3:4)'*Sol17.pdf_cond(3:4,:))')./sum(StaticParamsCC.p(3:4)),...
     (StaticParamsCC.kGrid*(StaticParamsCC.p(1:2)'*Sol18.pdf_cond(1:2,:))')./sum(StaticParamsCC.p(1:2)),...
     (StaticParamsCC.kGrid*(StaticParamsCC.p(3:4)'*Sol18.pdf_cond(3:4,:))')./sum(StaticParamsCC.p(3:4))];
Output(3,:) = [(StaticParams.kGrid*Sol1.pdf_cond'),...
     (StaticParams.kGrid*Sol5.pdf_cond'),...
     (StaticParams.kGrid*Sol6.pdf_cond'),...
     (StaticParams.kGrid*Sol7.pdf_cond'),...
     (StaticParams.kGrid*Sol11.pdf_cond'),...
     (StaticParams.kGrid*Sol12.pdf_cond'),...
     (StaticParamsCC.kGrid*Sol15.pdf_cond'),...
     (StaticParamsCC.kGrid*Sol17.pdf_cond'),...
     (StaticParamsCC.kGrid*Sol18.pdf_cond')];
     
% labor supply
Output(4,1:4:N) = [(StaticParams.lGrid*(pL'*Sol1.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol5.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol6.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol7.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol11.pdf_cond_l)'),...
     (StaticParams.lGrid*(pL'*Sol12.pdf_cond_l)'),...
     (StaticParamsCC.lGrid*(pL'*Sol15.pdf_cond_l)'),...
     (StaticParamsCC.lGrid*(pL'*Sol17.pdf_cond_l)'),...
     (StaticParamsCC.lGrid*(pL'*Sol18.pdf_cond_l)')];
Output(5,1:2:N) = [(StaticParams.lGrid*(pL(1:2)'*Sol1.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol1.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol5.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol5.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol6.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol6.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol7.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol7.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol11.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol11.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParams.lGrid*(pL(1:2)'*Sol12.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParams.lGrid*(pL(3:4)'*Sol12.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParamsCC.lGrid*(pL(1:2)'*Sol15.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParamsCC.lGrid*(pL(3:4)'*Sol15.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParamsCC.lGrid*(pL(1:2)'*Sol17.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParamsCC.lGrid*(pL(3:4)'*Sol17.pdf_cond_l(3:4,:))')./sum(pL(3:4)),...
     (StaticParamsCC.lGrid*(pL(1:2)'*Sol18.pdf_cond_l(1:2,:))')./sum(pL(1:2)),...
     (StaticParamsCC.lGrid*(pL(3:4)'*Sol18.pdf_cond_l(3:4,:))')./sum(pL(3:4))];  

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
Output(10,1:4:N) = [Sol1.agC*StaticParams.p,Sol5.agC*StaticParams.p,Sol6.agC*StaticParams.p,...
                    Sol7.agC*StaticParams.p,Sol11.agC*StaticParams.p,Sol12.agC*StaticParams.p,...
                    Sol15.agC*StaticParamsCC.p,Sol17.agC*StaticParamsCC.p,Sol18.agC*StaticParamsCC.p];
Output(11,1:2:N) = 2.*[Sol1.agC(1:2)*StaticParams.p(1:2),Sol1.agC(3:4)*StaticParams.p(3:4),...
                     Sol5.agC(1:2)*StaticParams.p(1:2),Sol5.agC(3:4)*StaticParams.p(3:4),...
                     Sol6.agC(1:2)*StaticParams.p(1:2),Sol6.agC(3:4)*StaticParams.p(3:4),...
                     Sol7.agC(1:2)*StaticParams.p(1:2),Sol7.agC(3:4)*StaticParams.p(3:4),...
                     Sol11.agC(1:2)*StaticParams.p(1:2),Sol11.agC(3:4)*StaticParams.p(3:4),...
                     Sol12.agC(1:2)*StaticParams.p(1:2),Sol12.agC(3:4)*StaticParams.p(3:4),...
                     Sol15.agC(1:2)*StaticParamsCC.p(1:2),Sol15.agC(3:4)*StaticParamsCC.p(3:4),...
                     Sol17.agC(1:2)*StaticParamsCC.p(1:2),Sol17.agC(3:4)*StaticParamsCC.p(3:4),...
                     Sol18.agC(1:2)*StaticParamsCC.p(1:2),Sol18.agC(3:4)*StaticParamsCC.p(3:4)];
Output(12,:) = [Sol1.agC,Sol5.agC,Sol6.agC,...
                Sol7.agC,Sol11.agC,Sol12.agC,...
                Sol15.agC,Sol17.agC,Sol18.agC];

% welfare
Output(13,1:4:N) = [Sol1.agWelfare*StaticParams.p,...
                    Sol5.agWelfare*StaticParams.p,...
                    Sol6.agWelfare*StaticParams.p,...
                    Sol7.agWelfare*StaticParams.p,...
                    Sol11.agWelfare*StaticParams.p,...
                    Sol12.agWelfare*StaticParams.p,...
                    Sol15.agWelfare*StaticParamsCC.p,...
                    Sol17.agWelfare*StaticParamsCC.p,...
                    Sol18.agWelfare*StaticParamsCC.p];
Output(14,1:2:N) = 2.*[Sol1.agWelfare(1:2)*StaticParams.p(1:2),Sol1.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol5.agWelfare(1:2)*StaticParams.p(1:2),Sol5.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol6.agWelfare(1:2)*StaticParams.p(1:2),Sol6.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol7.agWelfare(1:2)*StaticParams.p(1:2),Sol7.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol11.agWelfare(1:2)*StaticParams.p(1:2),Sol11.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol12.agWelfare(1:2)*StaticParams.p(1:2),Sol12.agWelfare(3:4)*StaticParams.p(3:4),...
                     Sol15.agWelfare(1:2)*StaticParamsCC.p(1:2),Sol15.agWelfare(3:4)*StaticParamsCC.p(3:4),...
                     Sol17.agWelfare(1:2)*StaticParamsCC.p(1:2),Sol17.agWelfare(3:4)*StaticParamsCC.p(3:4),...
                     Sol18.agWelfare(1:2)*StaticParamsCC.p(1:2),Sol18.agWelfare(3:4)*StaticParamsCC.p(3:4)];
Output(15,:) = [Sol1.agWelfare,Sol5.agWelfare,Sol6.agWelfare,...
                Sol7.agWelfare,Sol11.agWelfare,Sol12.agWelfare,...
                Sol15.agWelfare,Sol17.agWelfare,Sol18.agWelfare];
% derivatives
Output([1,4,6,8,10,13],6)=(Output([1,4,6,8,10,13],5)-Output([1,4,6,8,10,13],1))/norm(0.015*[1,1]);
Output([1,4,6,8,10,13],10)=(Output([1,4,6,8,10,13],9)-Output([1,4,6,8,10,13],1))/norm(0.015*[1,1]);
Output([1,4,6,8,10,13],18)=(Output([1,4,6,8,10,13],17)-Output([1,4,6,8,10,13],13))/norm(0.015*[1,0]);
Output([1,4,6,8,10,13],22)=(Output([1,4,6,8,10,13],21)-Output([1,4,6,8,10,13],13))/norm(0.015*[0,1]);
Output([1,4,6,8,10,13],30)=(Output([1,4,6,8,10,13],29)-Output([1,4,6,8,10,13],25))/norm(0.015*[1,0]);
Output([1,4,6,8,10,13],34)=(Output([1,4,6,8,10,13],33)-Output([1,4,6,8,10,13],25))/norm(0.015*[0,1]);

%--------------------------------------------------------------------------
% Fit welfare function to computed values, x grid: benfit in boom, 
% y grid: benefit in recession
%--------------------------------------------------------------------------
% case \bar{k}=0, 12 data points (folders 1,3,5:14)
%--------------------------------------------------------
fol = [1,3,5:14];
x = zeros(length(fol),1);
y = zeros(length(fol),1);
w = zeros(length(fol),1);
data_count = 0;
for case_nr=fol
    data_count = data_count+1;
    folder = strcat('res_case',num2str(case_nr));
    StaticParams = load(strcat(folder,'/StaticParams.mat'));
    switch case_nr
        case {ct{:}}
            count = 3;
        otherwise
            count = 2;
    end
    Sol = load(strcat(folder,'/Sol',num2str(count),'_PFI.mat'));
    x(data_count) = StaticParams.mu;
    y(data_count) = StaticParams.mu+StaticParams.mu_extra;
    w(data_count) = Sol.agWelfare*StaticParams.p;
end
[X,Y] = ndgrid(unique(x),unique(y));
Z = zeros(size(Y));
% fill with computed welfare values
for i=1:length(w)
    Z(logical((X==x(i)).*(Y==y(i))))=w(i);
end
% then interpolate along diagonal (acyclical policies)
ind = sub2ind(size(Z),1:size(Z,1),1:size(Z,2));
bool = Z(ind)>=0;
Z(ind(logical(bool))) = interp1(Y(ind(logical(1-bool))),Z(ind(logical(1-bool))),Y(ind(logical(bool))),'linear','extrap');
% then interpolate the rest
Z(1,Z(1,:)==0) = interp1(Y(1,Z(1,:)<0),Z(1,Z(1,:)<0),Y(1,Z(1,:)==0),'linear','extrap');
Z(Z(:,1)==0,1) = interp1(X(Z(:,1)<0,1),Z(Z(:,1)<0,1),X(Z(:,1)==0,1),'linear','extrap');
Z(end,Z(end,:)==0) = interp1(Y(end,Z(end,:)<0),Z(end,Z(end,:)<0),Y(end,Z(end,:)==0),'linear','extrap');
Z(Z(:,end)==0,end) = interp1(X(Z(:,end)<0,end),Z(Z(:,end)<0,end),X(Z(:,end)==0,end),'linear','extrap');
Z1 = Z;
Z2 = Z;
for i=2:size(Z,1)-1
    Z1(i,Z(i,:)==0) = interp1(Y(i,Z1(i,:)<0),Z1(i,Z1(i,:)<0),Y(i,Z1(i,:)==0),'linear','extrap');
end
for i=2:size(Z,2)-1
    Z2(Z2(:,i)==0,i) = interp1(X(Z2(:,i)<0,i),Z2(Z2(:,i)<0,i),X(Z2(:,i)==0,i),'linear','extrap');
end
Z = 0.5.*(Z1+Z2);
[X_large,Y_large] = ndgrid(linspace(0,0.15,101),linspace(0,0.15,101));
Z_large = interpn(unique(x),unique(y),Z,X_large,Y_large);


%--------------------------------------------------------
for case_nr = [15,17,18,35,37,38]
    folder = strcat('res_case',num2str(case_nr));
    StaticParams = load(strcat(folder,'/StaticParams.mat'));
    count = 2;
    Sol = load(strcat(folder,'/Sol',num2str(count),'_PFI.mat'));
    Poly = load(strcat(folder,'/Poly',num2str(count),'_PFI.mat'));
    weights0 = calcWeights(Sol.cdf_cond(1:2,:),Poly,StaticParams,0);
    weights1 = calcWeights(Sol.cdf_cond(3:4,:),Poly,StaticParams,1);
    c_prime0 = max(0,interp1(StaticParams.kGrid_pol,squeeze(...
                      [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,1,:)),Sol.idx,weights0);...
                       getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,2,:)),Sol.idx,weights0)])',StaticParams.kGrid)');
    c_prime0(c_prime0<eps^2) = 0;
    c_prime0(Sol.pdf_cond(1:2,:)==0) = max(eps,c_prime0(Sol.pdf_cond(1:2,:)==0));
    c_prime1 = max(0,interp1(StaticParams.kGrid_pol,squeeze(...
                      [getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,3,:)),Sol.idx,weights1);...
                       getCApprox(Sol.grid_N,Sol.distrGrid,squeeze(Sol.c(:,4,:)),Sol.idx,weights1)])',StaticParams.kGrid)');
    c_prime1(Sol.pdf_cond(3:4,:)==0) = max(eps,c_prime1(Sol.pdf_cond(3:4,:)==0));
    c_prime1(c_prime1<eps^2) = 0;
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
        ut0_b = log(max(eps^20,aggregator0));
        ut1_b = log(max(eps^20,aggregator1));
    else
        ut0 = (aggregator0).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
        ut1 = (aggregator1).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
        ut0_b = (max(eps^20,aggregator0)).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
        ut1_b = (max(eps^20,aggregator1)).^(1-StaticParams.gamma)./(1-StaticParams.gamma);
    end
    Sol.agC = [sum(c_prime0.*Sol.pdf_cond(1:2,:),2)',...
               sum(c_prime1.*Sol.pdf_cond(3:4,:),2)'];
    Sol.agWelfare = [sum(ut0.*Sol.pdf_cond(1:2,:),2)',...
                     sum(ut1.*Sol.pdf_cond(3:4,:),2)'];
    Sol.agWelfare_b = [sum(ut0_b.*Sol.pdf_cond(1:2,:),2)',...
                       sum(ut1_b.*Sol.pdf_cond(3:4,:),2)'];
    save(strcat(folder,'/Sol',num2str(count),'_PFI.mat'),'-struct','Sol'); 
end 

% case \bar{k}=-5, 10 data points (folders 15:24)
%--------------------------------------------------------
fol = 15:24;
x = zeros(length(fol),1);
y = zeros(length(fol),1);
w = zeros(length(fol),1);
w_b = zeros(length(fol),1);
data_count = 0;
for case_nr=fol
    data_count = data_count+1;
    folder = strcat('res_case',num2str(case_nr));
    StaticParams = load(strcat(folder,'/StaticParams.mat'));
    switch case_nr
        case {ct{:}}
            count = 3;
        otherwise
            count = 2;
    end
    Sol = load(strcat(folder,'/Sol',num2str(count),'_PFI.mat'));
    x(data_count) = StaticParams.mu;
    y(data_count) = StaticParams.mu+StaticParams.mu_extra;
    w(data_count) = Sol.agWelfare*StaticParams.p;
    w_b(data_count) = Sol.agWelfare_b*StaticParams.p;
end
[Xcc,Ycc] = ndgrid(unique(x),unique(y));
Zcc = zeros(size(Ycc));
Zcc_interp = zeros(size(Ycc));
% fill with computed welfare values
for i=1:length(w)
    Zcc(logical((Xcc==x(i)).*(Ycc==y(i))))=w(i);
    Zcc_interp(logical((Xcc==x(i)).*(Ycc==y(i))))=w_b(i);
end
% then interpolate along diagonal (acyclical policies)
ind = sub2ind(size(Zcc),1:size(Zcc,1),1:size(Zcc,2));
bool = Zcc(ind)>=0;
Zcc_interp(ind(logical(bool))) = interp1(Ycc(ind(logical(1-bool))),Zcc_interp(ind(logical(1-bool))),Ycc(ind(logical(bool))),'linear','extrap');
Zcc(ind(logical(bool))) = Zcc_interp(ind(logical(bool)));
% then interpolate the rest
Zcc(1,Zcc(1,:)==0) = interp1(Ycc(1,Zcc(1,:)<0),Zcc(1,Zcc(1,:)<0),Ycc(1,Zcc(1,:)==0),'linear','extrap');
Zcc(Zcc(:,1)==0,1) = interp1(Xcc(Zcc(:,1)<0,1),Zcc(Zcc(:,1)<0,1),Xcc(Zcc(:,1)==0,1),'linear','extrap');
Zcc_interp(1,Zcc_interp(1,:)==0) = interp1(Ycc(1,Zcc_interp(1,:)<0),Zcc_interp(1,Zcc_interp(1,:)<0),Ycc(1,Zcc_interp(1,:)==0),'linear','extrap');
Zcc_interp(Zcc_interp(:,1)==0,1) = interp1(Xcc(Zcc_interp(:,1)<0,1),Zcc_interp(Zcc_interp(:,1)<0,1),Xcc(Zcc_interp(:,1)==0,1),'linear','extrap');
Zcc_interp(end,Zcc_interp(end,:)==0) = interp1(Ycc(end,Zcc_interp(end,:)<0),Zcc_interp(end,Zcc_interp(end,:)<0),Ycc(end,Zcc_interp(end,:)==0),'linear','extrap');
Zcc(end,Zcc(end,:)==0) = Zcc_interp(end,Zcc(end,:)==0);
Zcc_interp(Zcc_interp(:,end)==0,end) = interp1(Xcc(Zcc_interp(:,end)<0,end),Zcc_interp(Zcc_interp(:,end)<0,end),Xcc(Zcc_interp(:,end)==0,end),'linear','extrap');
Zcc(Zcc(:,end)==0,end) = Zcc_interp(Zcc(:,end)==0,end);
Z1 = Zcc_interp;
Z2 = Zcc_interp;
for i=2:size(Zcc,1)-1
    Z1(i,Z1(i,:)==0) = interp1(Ycc(i,Zcc_interp(i,:)<0),Zcc_interp(i,Zcc_interp(i,:)<0),Ycc(i,Zcc_interp(i,:)==0),'linear','extrap');
end
for i=2:size(Zcc,2)-1
    Z2(Z2(:,i)==0,i) = interp1(Xcc(Zcc_interp(:,i)<0,i),Zcc_interp(Zcc_interp(:,i)<0,i),Xcc(Zcc_interp(:,i)==0,i),'linear','extrap');
end
Zcc_interp = 0.5.*(Z1+Z2);
Zcc(Zcc==0) = Zcc_interp(Zcc==0);
[Xcc_large,Ycc_large] = ndgrid(linspace(0,0.15,101),linspace(0,0.15,101));
Zcc_large = interpn(unique(x),unique(y),max(log(eps^20),Zcc),Xcc_large,Ycc_large,'linear');

% case \bar{k}=-5(recession)/-4.5(boom), 10 data points (folders 35:44)
%--------------------------------------------------------
fol = 35:44;
x = zeros(length(fol),1);
y = zeros(length(fol),1);
w = zeros(length(fol),1);
w_b = zeros(length(fol),1);
data_count = 0;
for case_nr=fol
    data_count = data_count+1;
    folder = strcat('res_case',num2str(case_nr));
    StaticParams = load(strcat(folder,'/StaticParams.mat'));
    switch case_nr
        case {ct{:}}
            count = 3;
        otherwise
            count = 2;
    end
    Sol = load(strcat(folder,'/Sol',num2str(count),'_PFI.mat'));
    x(data_count) = StaticParams.mu;
    y(data_count) = StaticParams.mu+StaticParams.mu_extra;
    w(data_count) = Sol.agWelfare*StaticParams.p;
    w_b(data_count) = Sol.agWelfare_b*StaticParams.p;
end
[Xcc2,Ycc2] = ndgrid(unique(x),unique(y));
Zcc2 = zeros(size(Ycc2));
Zcc2_interp = zeros(size(Ycc2));
% fill with computed welfare values
for i=1:length(w)
    Zcc2(logical((Xcc2==x(i)).*(Ycc2==y(i))))=w(i);
    Zcc2_interp(logical((Xcc2==x(i)).*(Ycc2==y(i))))=w_b(i);
end
% then interpolate along diagonal (acyclical policies)
ind = sub2ind(size(Zcc2),1:size(Zcc2,1),1:size(Zcc2,2));
bool = Zcc2(ind)>=0;
Zcc2_interp(ind(logical(bool))) = interp1(Ycc2(ind(logical(1-bool))),Zcc2_interp(ind(logical(1-bool))),Ycc2(ind(logical(bool))),'linear','extrap');
Zcc2(ind(logical(bool))) = Zcc2_interp(ind(logical(bool)));
% then interpolate the rest
Zcc2_interp(1,Zcc2_interp(1,:)==0) = interp1(Ycc2(1,Zcc2_interp(1,:)<0),Zcc2_interp(1,Zcc2_interp(1,:)<0),Ycc2(1,Zcc2_interp(1,:)==0),'linear','extrap');
Zcc2_interp(Zcc2_interp(:,1)==0,1) = interp1(Xcc2(Zcc2_interp(:,1)<0,1),Zcc2_interp(Zcc2_interp(:,1)<0,1),Xcc2(Zcc2_interp(:,1)==0,1),'linear','extrap');
Zcc2(1,Zcc2(1,:)==0) = Zcc2_interp(1,Zcc2(1,:)==0);
Zcc2(Zcc2(:,1)==0,1) = interp1(Xcc2(Zcc2(:,1)<0,1),Zcc2(Zcc2(:,1)<0,1),Xcc2(Zcc2(:,1)==0,1),'linear','extrap');
Zcc2_interp(end,Zcc2_interp(end,:)==0) = interp1(Ycc2(end,Zcc2_interp(end,:)<0),Zcc2_interp(end,Zcc2_interp(end,:)<0),Ycc2(end,Zcc2_interp(end,:)==0),'linear','extrap');
Zcc2(end,Zcc2(end,:)==0) = Zcc2_interp(end,Zcc2(end,:)==0);
Zcc2_interp(Zcc2_interp(:,end)==0,end) = interp1(Xcc2(Zcc2_interp(:,end)<0,end),Zcc2_interp(Zcc2_interp(:,end)<0,end),Xcc2(Zcc2_interp(:,end)==0,end),'linear','extrap');
Zcc2(Zcc2(:,end)==0,end) = Zcc2_interp(Zcc2(:,end)==0,end);
Z1 = Zcc2_interp;
Z2 = Zcc2_interp;
for i=2:size(Zcc2,1)-1
    Z1(i,Z1(i,:)==0) = interp1(Ycc2(i,Zcc2_interp(i,:)<0),Zcc2_interp(i,Zcc2_interp(i,:)<0),Ycc2(i,Zcc2_interp(i,:)==0),'linear','extrap');
end
for i=2:size(Zcc2,2)-1
    Z2(Z2(:,i)==0,i) = interp1(Xcc2(Zcc2_interp(:,i)<0,i),Zcc2_interp(Zcc2_interp(:,i)<0,i),Xcc2(Zcc2_interp(:,i)==0,i),'linear','extrap');
end
Zcc2_interp = 0.5.*(Z1+Z2);
Zcc2(Zcc2==0) = Zcc2_interp(Zcc2==0);
Zcc2_large = interpn(unique(x),unique(y),max(log(eps^20),Zcc2),Xcc_large,Ycc_large,'linear');

%--------------------------------------------------------
figure(6)
subplot(1,3,1)
hold on
v=linspace(min(min(Z_large)),max(max(Z_large)),11);
C=contour(X_large,Y_large,Z_large,v);
colormap(winter)
idx_x=1:size(C,2);
idx_x=idx_x((C(1,:)==0.015));
idx_y=unique(C(2,idx_x));
text(0.015*ones(1,7),idx_y,string(round(v(10:-1:4),4)))
plot(X_large(Z_large==max(max(Z_large))),Y_large(Z_large==max(max(Z_large))),'r*')
hold off
axis([0 0.15 0 0.15])
xlabel('\nu in boom');
ylabel('\nu in recession');
title('Panel A')

subplot(1,3,2)
hold on
v=unique([linspace(min(min(Zcc_large)),min(min(Zcc_interp)),4),linspace(min(min(Zcc_interp)),max(max(Zcc_interp)),11)]);
C=contour(Xcc_large,Ycc_large,Zcc_large,v);
colormap(winter)
idx_x=1:size(C,2);
idx_x=idx_x((C(1,:)==0.015));
idx_y=unique(C(2,idx_x));
text(0.015*ones(1,10),idx_y([1,13:21]),string(round(v([1,13:-1:5]),4)))
plot(Xcc_large(Zcc_large==max(max(Zcc_large))),Ycc_large(Zcc_large==max(max(Zcc_large))),'r*')
hold off
xlabel('\nu in boom');
ylabel('\nu in recession');
title('Panel B')

subplot(1,3,3)
hold on
v=unique([linspace(min(min(Zcc2_large)),min(min(Zcc2_interp)),4),linspace(min(min(Zcc2_interp)),max(max(Zcc2_interp)),11)]);
C=contour(Xcc_large,Ycc_large,Zcc2_large,v);
colormap(winter)
idx_x=1:size(C,2);
idx_x=idx_x((C(1,:)==0.015));
idx_y=unique(C(2,idx_x));
text(0.015*ones(1,6),idx_y([1,13:17]),string(round(v([1,13:-1:9]),4)))
plot(Xcc_large(Zcc2_large==max(max(Zcc2_large))),Ycc_large(Zcc2_large==max(max(Zcc2_large))),'r*')
hold off
xlabel('\nu in boom');
ylabel('\nu in recession');
title('Panel C')