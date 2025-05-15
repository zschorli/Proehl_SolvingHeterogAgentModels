% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file produces the plots in the paper.
%__________________________________________________________________________
clc;
clear all;

case_nr = 1;
folder = strcat('res_case',num2str(case_nr));
StaticParams = load(strcat(folder,'/StaticParams.mat'));

%--------------------------------------------------------------------------
% Base distribution
%--------------------------------------------------------------------------
Poly =load(strcat(folder,'/Poly3_PPA.mat'));
figure(1)
P=Poly.pdf{2};
MP=[1,(P(1,2:end-2)>P(1,3:end-1)).*(P(1,2:end-2)>P(1,1:end-3)),0];
MP(1,P(1,:)==max(P(1,:)))=0;
kGrid_midPts = [StaticParams.kGrid(1),(StaticParams.kGrid(2:end)+StaticParams.kGrid(1:end-1))/2];
subplot(1,2,1)
hold on
plot(sort([StaticParams.kGrid,StaticParams.kGrid(1:end-1)]),P(1,sort([1:length(P),2:length(P)])),'b');
plot(kGrid_midPts(logical(MP(1,1:end))),P(1,logical(MP(1,:))),'rs');
axis([0 140 0 3e-3])
legend('histogram of the density','mass points');
xlabel('individual capital k');
ylabel('probability');
title('Panel A')
subplot(1,2,2)
hold on
plot(sort([StaticParams.kGrid,StaticParams.kGrid(1:end-1)]),P(1,sort([1:length(P),2:length(P)])),'b');
plot(StaticParams.kGrid([1 1]),P(1,1:2),'b','LineWidth',1);
plot(kGrid_midPts(logical(MP(1,1:end))),P(1,logical(MP(1,:))),'rs');
axis([0 5 0 2e-4])
xlabel('individual capital k');
ylabel('probability');
title('Panel B')

%--------------------------------------------------------------------------
% Orthogonal Polynomials
%--------------------------------------------------------------------------
figure(2)
subplot(1,3,1)
hold on
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{2}(1,:)*Poly.xi{2},'b')
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{2}(2,:)*Poly.xi{2},'g')
axis([0 60 -40 20])
legend('\Phi^k_0','\Phi^k_1','Location','SouthEast');
xlabel('basic random variable \xi^k');
ylabel('function value of the polynomial');
title('Panel A')
subplot(1,3,2)
hold on
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{2}(3,:)*Poly.xi{2},'r')
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{2}(4,:)*Poly.xi{2},'c')
axis([0 150 -2.5e4 2.5e4])
legend('\Phi^k_2','\Phi^k_3','Location','SouthEast');
xlabel('basic random variable \xi^k');
ylabel('function value of the polynomial');
title('Panel B')
subplot(1,3,3)
hold on
plot(Poly.xi{1}(2,:),Poly.phi_coeffs{1}(1,:)*Poly.xi{1},'b*')
plot(Poly.xi{1}(2,:),Poly.phi_coeffs{1}(2,:)*Poly.xi{1},'g*')
plot(Poly.xi{1}(2,:),Poly.phi_coeffs{1}(3,:)*Poly.xi{1},'r*')
plot(Poly.xi{1}(2,:),Poly.phi_coeffs{1}(4,:)*Poly.xi{1},'c*')
hold off
axis([0.5 3.5 -2 1.2])
legend('\Phi^z_0','\Phi^z_1','\Phi^z_2','\Phi^z_3','Location','SouthEast');
xlabel('basic random variable \xi^z');
ylabel('function value of the polynomial');
title('Panel C')


%--------------------------------------------------------------------------
% Interpretation of Polynomials
%--------------------------------------------------------------------------
figure(3)
hold on
P=Poly.pdf{2};
Pol = @(n,k) interp1(Poly.xi{2}(2,:),Poly.phi_coeffs{2}(n+1,:)*Poly.xi{2},k);
[a1,b1] = sort(36+1*Pol(1,StaticParams.kGrid));
a1 = interp1(a1,cumsum(P(1,b1)),max(a1(1),StaticParams.kGrid));
b1 = [a1(1),a1(2:end)-a1(1:end-1)];
[a2,b2] = sort(36+1*Pol(1,StaticParams.kGrid)+5e-3*Pol(2,StaticParams.kGrid));
a2 = interp1(a2,cumsum(P(1,b2)),max(a2(1),StaticParams.kGrid));
b2 = [a2(1),a2(2:end)-a2(1:end-1)];
[a3,b3] = sort(36+1*Pol(1,StaticParams.kGrid)+5e-3*Pol(2,StaticParams.kGrid)+5e-5*Pol(3,StaticParams.kGrid));
a3 = interp1(a3,cumsum(P(1,b3)),max(a3(1),StaticParams.kGrid));
b3 = [a3(1),a3(2:end)-a3(1:end-1)];
g = sort([StaticParams.kGrid,StaticParams.kGrid(2:end)-1e-15]);
plot([36-1e-16 36 36+1e-16],[0,1,0],'b','LineWidth',1.5)
plot(g,b1(sort([1:length(b1),1:length(b1)-1])),'g')
plot(g,b2(sort([1:length(b1),1:length(b1)-1])),'r')
plot(g,b3(sort([1:length(b1),1:length(b1)-1])),'c')
legend('36*\Phi^k_0','...+1*\Phi^k_1','...+5e-3*\Phi^k_2','...+5e-5*\Phi^k_3');
axis([0 80 0 4e-3])
xlabel('individual capital k');
ylabel('probability');

for case_nr = 1:3
switch case_nr
    case 1
        folder = 'res_case1';
    case 2
        folder = 'res_case2';
    case 3
        folder = 'res_case3';
end
StaticParams = load(strcat(folder,'/StaticParams.mat'));

%--------------------------------------------------------------------------
% Boxplot of Euler equation errors
%--------------------------------------------------------------------------
Res = load(strcat(folder,'/Res.mat'));
a1=[Res.stErr_KS1;Res.stErr_KS2;Res.stErr_KS3;Res.stErr_KS4;...
   Res.stErr_PPA0;Res.stErr_PPA1;Res.stErr_PPA2;Res.stErr_PPA3;...
   Res.stErr_PFI0;Res.stErr_PFI1;Res.stErr_PFI2;Res.stErr_PFI3];
b1=[zeros(length(Res.stErr_KS1),1);1*ones(length(Res.stErr_KS2),1);...
   2*ones(length(Res.stErr_KS3),1);3*ones(length(Res.stErr_KS4),1);...
   4*ones(length(Res.stErr_PPA0),1);5*ones(length(Res.stErr_PPA1),1);...
   6*ones(length(Res.stErr_PPA2),1);7*ones(length(Res.stErr_PPA3),1);...
   ...8*ones(length(Res.stErr_PPA4),1);
   9*ones(length(Res.stErr_PFI0),1);...
   10*ones(length(Res.stErr_PFI1),1);11*ones(length(Res.stErr_PFI2),1);...
   12*ones(length(Res.stErr_PFI3),1)];%13*ones(length(Res.stErr_PFI4),1)];
a2=[Res.dynErr_KS1;Res.dynErr_KS2;Res.dynErr_KS3;Res.dynErr_KS4;...
   Res.dynErr_PPA0;Res.dynErr_PPA1;Res.dynErr_PPA2;Res.dynErr_PPA3;...
   Res.dynErr_PFI0;Res.dynErr_PFI1;Res.dynErr_PFI2;Res.dynErr_PFI3];
b2=[zeros(length(Res.dynErr_KS1),1);1*ones(length(Res.dynErr_KS2),1);...
   2*ones(length(Res.dynErr_KS3),1);3*ones(length(Res.dynErr_KS4),1);...
   4*ones(length(Res.dynErr_PPA0),1);5*ones(length(Res.dynErr_PPA1),1);...
   6*ones(length(Res.dynErr_PPA2),1);7*ones(length(Res.dynErr_PPA3),1);...
   ...8*ones(length(Res.dynErr_PPA4),1);
   9*ones(length(Res.dynErr_PFI0),1);...
   10*ones(length(Res.dynErr_PFI1),1);11*ones(length(Res.dynErr_PFI2),1);...
   12*ones(length(Res.dynErr_PFI3),1)];%;13*ones(length(Res.dynErr_PFI4),1)];
figure(4*(case_nr-1)+4)
subplot(1,2,1)
boxplot(a1,b1,...
         'Labels',{'K-S1','K-S2','K-S3','K-S4','PPA0','PPA1','PPA2','PPA3','PFI0','PFI1','PFI2','PFI3'},...
         'OutlierSize',1,'Symbol','r.');%'Whisker',4
title('Panel A')
axis([0.5 12.5 -9 -1])
subplot(1,2,2)
boxplot(a2,b2,...
         'Labels',{'K-S1','K-S2','K-S3','K-S4','PPA0','PPA1','PPA2','PPA3','PFI0','PFI1','PFI2','PFI3'},...
         'OutlierSize',1,'Symbol','r.');%'Whisker',4
title('Panel B')   
axis([0.5 12.5 -9 -1])
clear a1 b1 a2 b2;


%--------------------------------------------------------------------------
% Boxplot of law of motion errors
%--------------------------------------------------------------------------
a1=[Res.LoMErr_KS1;Res.LoMErr_KS2;Res.LoMErr_KS3;Res.LoMErr_KS4;...
    Res.LoMErr_PPA0;Res.LoMErr_PPA1;Res.LoMErr_PPA2;Res.LoMErr_PPA3;...Res.LoMErr_PPA4;...
    Res.LoMErr_PFI0;Res.LoMErr_PFI1;Res.LoMErr_PFI2;Res.LoMErr_PFI3];%Res.LoMErr_PFI4];
b1=[zeros(length(Res.LoMErr_KS1),1);ones(length(Res.LoMErr_KS2),1);...
   2*ones(length(Res.LoMErr_KS3),1);3*ones(length(Res.LoMErr_KS4),1);...
   4*ones(length(Res.LoMErr_PPA0),1);5*ones(length(Res.LoMErr_PPA1),1);...
   6*ones(length(Res.LoMErr_PPA2),1);7*ones(length(Res.LoMErr_PPA3),1);...
   ...8*ones(length(Res.LoMErr_PPA4),1);
   9*ones(length(Res.LoMErr_PFI0),1);...
   10*ones(length(Res.LoMErr_PFI1),1);11*ones(length(Res.LoMErr_PFI2),1);...
   12*ones(length(Res.LoMErr_PFI3),1)];%13*ones(length(Res.LoMErr_PFI4),1)];
figure(4*(case_nr-1)+5)
boxplot(a1,b1,...
         'Labels',{'K-S1','K-S2','K-S3','K-S4','PPA0','PPA1','PPA2','PPA3','PFI0','PFI1','PFI2','PFI3'},...
         'OutlierSize',1,'Symbol','r.');%'Whisker',4
axis([0.5 12.5 -8 -2])


%--------------------------------------------------------------------------
% Stationary distributions 
%--------------------------------------------------------------------------
Sol0 = load(strcat(folder,'/Sol0_PPA.mat'));
Sol1 = load(strcat(folder,'/Sol1_PPA.mat'));
Sol2 = load(strcat(folder,'/Sol2_PPA.mat'));
Sol3 = load(strcat(folder,'/Sol3_PPA.mat'));
% Sol4 = load(strcat(folder,'/Sol4_PFI.mat'));
KS = load(strcat(folder,'/KS1_Sol',num2str(case_nr),'.mat'));
KS2 = load(strcat(folder,'/KS2_Sol',num2str(case_nr),'.mat'));
KS3 = load(strcat(folder,'/KS3_Sol',num2str(case_nr),'.mat'));
KS4 = load(strcat(folder,'/KS4_Sol',num2str(case_nr),'.mat'));

p=[sum(StaticParams.p([1,3])),sum(StaticParams.p([2,4]))];

figure(4*(case_nr-1)+6)
hold on
plot(StaticParams.kGrid,p*(Sol1.pdf-Sol0.pdf),'b')
plot(StaticParams.kGrid,p*(Sol2.pdf-Sol1.pdf),'r')
plot(StaticParams.kGrid,p*(Sol3.pdf-Sol2.pdf),'g')
% plot(StaticParams.kGrid,p*(Sol4.pdf-Sol3.pdf),'m')
plot([0 150],StaticParams.criter_k*ones(1,2),'k',[0 150],-StaticParams.criter_k*ones(1,2),'k')
hold off 
axis([0 150 -3e-4 3e-4])
legend('1^{st} order - 0^{th} order','2^{nd} order - 1^{st} order','3^{rd} order - 2^{nd} order', 'termination criterion','Location','northeast')
xlabel('individual capital k');
ylabel('change in probability');

figure(4*(case_nr-1)+7)
subplot(1,2,1)
hold on
plot(StaticParams.kGrid,p*KS.pdf,'k','LineWidth',1);
plot(StaticParams.kGrid,p*Sol0.pdf,'r','LineWidth',1);
plot(StaticParams.kGrid,p*Sol1.pdf,'c','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol3.pdf,'m--','LineWidth',1.5);
% plot(StaticParams.kGrid,p*Sol4.pdf,'g-.','LineWidth',1.5);
hold off
axis([0 100 0 0.0035])
% axis([0 100 0 0.012])
legend('K-S 1 moment','PPA 0^{th} order','PPA 1^{st} order','PPA 2^{nd} order','PPA 3^{rd} order')
xlabel('individual capital k');
ylabel('probability');
title('Panel A')
subplot(1,2,2)
hold on
plot(StaticParams.kGrid,p*Sol3.pdf,'k','LineWidth',1);
plot(StaticParams.kGrid,p*KS.pdf,'r','LineWidth',1);
plot(StaticParams.kGrid,p*KS2.pdf,'c','LineWidth',1.5);
plot(StaticParams.kGrid,p*KS3.pdf,'b','LineWidth',1.5);
plot(StaticParams.kGrid,p*KS4.pdf,'m--','LineWidth',1.5);
hold off
axis([0 100 0 0.0035])
% axis([0 100 0 0.012])
legend('PPA 3^{rd} order','K-S 1 moment','K-S 2 moments','K-S 3 moments','K-S 4 moments')
xlabel('individual capital k');
ylabel('probability');
title('Panel B')
end

%--------------------------------------------------------------------------
% High vs. low unempl. benefit
%--------------------------------------------------------------------------

StaticParams1 = load(strcat('res_case1/StaticParams.mat'));
StaticParams2 = load(strcat('res_case2/StaticParams.mat'));
StaticParams3 = load(strcat('res_case3/StaticParams.mat'));
preSol1 = load(strcat('res_case1/preSol_PPA.mat'));
preSol2 = load(strcat('res_case2/preSol_PPA.mat'));
preSol3 = load(strcat('res_case3/preSol_PPA.mat'));
Sol1 = load(strcat('res_case1/Sol3_PPA.mat'));
Sol2 = load(strcat('res_case2/Sol3_PPA.mat'));
Sol3 = load(strcat('res_case3/Sol3_PPA.mat'));
compSol = load(strcat('compareDistr_NoAggShock/compareSol.mat'));

figure(13)
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
axis([0 100 0 0.03])
legend('\nu=15%','\nu=40%','\nu=65%')
xlabel('individual capital k');
ylabel('probability');

figure(14)
pdf = compSol.pdf_k;
pdf(:,2:end) = pdf(:,2:end)./repmat(StaticParams1.kGrid(2:end)-StaticParams1.kGrid(1:end-1),[size(pdf,1),1]); 
plot(StaticParams1.kGrid,pdf);
axis([0 100 0 0.03])
names = cell(1,length(compSol.muGrid));
for i=1:length(compSol.muGrid)
    names{i} = strcat('\nu=',num2str(compSol.muGrid(i)*100),'%');
end
legend(names)
xlabel('individual capital k');
ylabel('probability');


figure(15)
subplot(1,3,1)
hold on 
pdf = 2*StaticParams1.p(1:2)'*preSol1.pdf_cond_b;
pdf(2:end) = pdf(2:end)./(StaticParams1.kGrid(2:end)-StaticParams1.kGrid(1:end-1)); 
plot(StaticParams1.kGrid,pdf,'r')
pdf = 2*StaticParams1.p(3:4)'*preSol1.pdf_cond_g;
pdf(2:end) = pdf(2:end)./(StaticParams1.kGrid(2:end)-StaticParams1.kGrid(1:end-1)); 
plot(StaticParams1.kGrid,pdf,'b')
pdf = p*Sol1.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams1.kGrid(2:end)-StaticParams1.kGrid(1:end-1)); 
plot(StaticParams1.kGrid,pdf,'g','LineWidth',1.5);
hold off
axis([0 140 0 0.026])
legend('no aggr. risk: bad state','no aggr. risk: good state','aggr. risk with both states')
xlabel('individual capital k');
ylabel('probability');
title('Panel A')%\nu=15%

subplot(1,3,2)
hold on 
pdf = 2*StaticParams2.p(1:2)'*preSol2.pdf_cond_b;
pdf(2:end) = pdf(2:end)./(StaticParams2.kGrid(2:end)-StaticParams2.kGrid(1:end-1)); 
plot(StaticParams2.kGrid,pdf,'r')
pdf = 2*StaticParams2.p(3:4)'*preSol2.pdf_cond_g;
pdf(2:end) = pdf(2:end)./(StaticParams2.kGrid(2:end)-StaticParams2.kGrid(1:end-1)); 
plot(StaticParams2.kGrid,pdf,'b')
pdf = p*Sol2.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams2.kGrid(2:end)-StaticParams2.kGrid(1:end-1)); 
plot(StaticParams2.kGrid,pdf,'g','LineWidth',1.5);
hold off
axis([0 140 0 0.026])
xlabel('individual capital k');
ylabel('probability');
title('Panel B')%\nu=40%

subplot(1,3,3)
hold on 
pdf = 2*StaticParams3.p(1:2)'*preSol3.pdf_cond_b;
pdf(2:end) = pdf(2:end)./(StaticParams3.kGrid(2:end)-StaticParams3.kGrid(1:end-1)); 
plot(StaticParams3.kGrid,pdf,'r')
pdf = 2*StaticParams3.p(3:4)'*preSol3.pdf_cond_g;
pdf(2:end) = pdf(2:end)./(StaticParams3.kGrid(2:end)-StaticParams3.kGrid(1:end-1)); 
plot(StaticParams3.kGrid,pdf,'b')
pdf = p*Sol3.pdf;
pdf(2:end) = pdf(2:end)./(StaticParams3.kGrid(2:end)-StaticParams3.kGrid(1:end-1)); 
plot(StaticParams3.kGrid,pdf,'g','LineWidth',1.5);
hold off
axis([0 140 0 0.026])
xlabel('individual capital k');
ylabel('probability');
title('Panel C')%\nu=65%