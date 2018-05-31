% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
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
axis([0 140 0 2.5e-3])
legend('histogram of the density','mass points');
xlabel('individual capital k');
ylabel('probability');
title('Panel A')
subplot(1,2,2)
hold on
plot(sort([StaticParams.kGrid,StaticParams.kGrid(1:end-1)]),P(1,sort([1:length(P),2:length(P)])),'b');
plot(StaticParams.kGrid([1 1]),P(1,1:2),'b','LineWidth',1);
plot(kGrid_midPts(logical(MP(1,1:end))),P(1,logical(MP(1,:))),'rs');
axis([0 5 0 1.2e-4])
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
axis([0 120 -1.6e4 1e4])
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
[a2,b2] = sort(36+1*Pol(1,StaticParams.kGrid)+1e-2*Pol(2,StaticParams.kGrid));
a2 = interp1(a2,cumsum(P(1,b2)),max(a2(1),StaticParams.kGrid));
b2 = [a2(1),a2(2:end)-a2(1:end-1)];
[a3,b3] = sort(36+1*Pol(1,StaticParams.kGrid)+1e-2*Pol(2,StaticParams.kGrid)+1e-4*Pol(3,StaticParams.kGrid));
a3 = interp1(a3,cumsum(P(1,b3)),max(a3(1),StaticParams.kGrid));
b3 = [a3(1),a3(2:end)-a3(1:end-1)];
g = sort([StaticParams.kGrid,StaticParams.kGrid(2:end)-1e-15]);
plot([36-1e-16 36 36+1e-16],[0,1,0],'b','LineWidth',1.5)
plot(g,b1(sort([1:length(b1),1:length(b1)-1])),'g')
plot(g,b2(sort([1:length(b1),1:length(b1)-1])),'r')
plot(g,b3(sort([1:length(b1),1:length(b1)-1])),'c')
legend('36*\Phi^k_0','...+1*\Phi^k_1','...+1e-2*\Phi^k_2','...+1e-4*\Phi^k_3');
axis([0 80 0 5e-3])
xlabel('individual capital k');
ylabel('probability');

for case_nr = 1:2
switch case_nr
    case 1
        folder = 'res_case1';
    case 2
        folder = 'res_case2';
end
StaticParams = load(strcat(folder,'/StaticParams.mat'));
%--------------------------------------------------------------------------
% Boxplot of Euler equation errors
%--------------------------------------------------------------------------
Res = load(strcat(folder,'/Res.mat'));
if case_nr==1
    a1=[Res.stErr_KS1;Res.stErr_R;Res.stErr_DR;...
       Res.stErr_PPA0;Res.stErr_PPA1;Res.stErr_PPA2;Res.stErr_PPA3;...
       Res.stErr_PFI0;Res.stErr_PFI1;Res.stErr_PFI2;Res.stErr_PFI3];
    b1=[zeros(length(Res.stErr_KS1),1);ones(length(Res.stErr_R),1);...
       2*ones(length(Res.stErr_DR),1);3*ones(length(Res.stErr_PPA0),1);...
       4*ones(length(Res.stErr_PPA1),1);5*ones(length(Res.stErr_PPA2),1);...
       6*ones(length(Res.stErr_PPA3),1);7*ones(length(Res.stErr_PFI0),1);...
       8*ones(length(Res.stErr_PFI1),1);9*ones(length(Res.stErr_PFI2),1);...
       10*ones(length(Res.stErr_PFI3),1)];
    a2=[Res.dynErr_KS1;Res.dynErr_R;Res.dynErr_DR;...
       Res.dynErr_PPA0;Res.dynErr_PPA1;Res.dynErr_PPA2;Res.dynErr_PPA3;...
       Res.dynErr_PFI0;Res.dynErr_PFI1;Res.dynErr_PFI2;Res.dynErr_PFI3];
    b2=[zeros(length(Res.dynErr_KS1),1);ones(length(Res.dynErr_R),1);...
       2*ones(length(Res.dynErr_DR),1);3*ones(length(Res.dynErr_PPA0),1);...
       4*ones(length(Res.dynErr_PPA1),1);5*ones(length(Res.dynErr_PPA2),1);...
       6*ones(length(Res.dynErr_PPA3),1);7*ones(length(Res.dynErr_PFI0),1);...
       8*ones(length(Res.dynErr_PFI1),1);9*ones(length(Res.dynErr_PFI2),1);...
       10*ones(length(Res.dynErr_PFI3),1)];  
else
    a1=[Res.stErr_KS1;Res.stErr_R;Res.stErr_DR;...
       Res.stErr_PPA0;Res.stErr_PPA2;Res.stErr_PPA3;Res.stErr_PPA4;...
       Res.stErr_PFI0;Res.stErr_PFI1;Res.stErr_PFI2;Res.stErr_PFI3;Res.stErr_PFI4];
    b1=[zeros(length(Res.stErr_KS1),1);ones(length(Res.stErr_R),1);...
       2*ones(length(Res.stErr_DR),1);3*ones(length(Res.stErr_PPA0),1);...
       4*ones(length(Res.stErr_PPA2),1);5*ones(length(Res.stErr_PPA3),1);...
       6*ones(length(Res.stErr_PPA4),1);7*ones(length(Res.stErr_PFI0),1);...
       8*ones(length(Res.stErr_PFI1),1);9*ones(length(Res.stErr_PFI2),1);...
       10*ones(length(Res.stErr_PFI3),1);11*ones(length(Res.stErr_PFI4),1)];
    a2=[Res.dynErr_KS1;Res.dynErr_R;Res.dynErr_DR;...
       Res.dynErr_PPA0;Res.dynErr_PPA2;Res.dynErr_PPA3;Res.dynErr_PPA4;...
       Res.dynErr_PFI0;Res.dynErr_PFI1;Res.dynErr_PFI2;Res.dynErr_PFI3;Res.dynErr_PFI4];
    b2=[zeros(length(Res.dynErr_KS1),1);ones(length(Res.dynErr_R),1);...
       2*ones(length(Res.dynErr_DR),1);3*ones(length(Res.dynErr_PPA0),1);...
       4*ones(length(Res.dynErr_PPA2),1);5*ones(length(Res.dynErr_PPA3),1);...
       6*ones(length(Res.dynErr_PPA4),1);7*ones(length(Res.dynErr_PFI0),1);...
       8*ones(length(Res.dynErr_PFI1),1);9*ones(length(Res.dynErr_PFI2),1);...
       10*ones(length(Res.dynErr_PFI3),1);11*ones(length(Res.dynErr_PFI4),1)];
end
figure(4*(case_nr-1)+4)
subplot(1,2,1)
if case_nr==1
    boxplot(abs(a1),b1,...
             'Labels',{'K-S','R','D-R','PPA0','PPA1','PPA2','PPA3','PFI0','PFI1','PFI2','PFI3'},...
             'Whisker',4,'OutlierSize',1,'Symbol','r.');
else
    boxplot(abs(a1),b1,...
         'Labels',{'K-S','R','D-R','PPA0','PPA2','PPA3','PPA4','PFI0','PFI1','PFI2','PFI3','PFI4'},...
         'Whisker',4,'OutlierSize',1,'Symbol','r.');
end
title('Panel A')
if case_nr==1
    axis([0.5 11.5 -5e-5 9e-4])
else
    axis([0.5 12.5 -5e-5 9e-4])
end
subplot(1,2,2)
if case_nr==1
    boxplot(abs(a2),b2,...
             'Labels',{'K-S','R','D-R','PPA0','PPA1','PPA2','PPA3','PFI0','PFI1','PFI2','PFI3'},...
             'Whisker',4,'OutlierSize',1,'Symbol','r.');
else
    boxplot(abs(a2),b2,...
         'Labels',{'K-S','R','D-R','PPA0','PPA2','PPA3','PPA4','PFI0','PFI1','PFI2','PFI3','PFI4'},...
         'Whisker',4,'OutlierSize',1,'Symbol','r.');
end
title('Panel B')   
if case_nr==1
    axis([0.5 11.5 -5e-5 9e-4])
else
    axis([0.5 12.5 -5e-5 9e-4])
end

figure(4*(case_nr-1)+5)
subplot(1,2,1)
if case_nr==1
    boxplot(a1,b1,...
             'Labels',{'K-S','R','D-R','PPA0','PPA1','PPA2','PPA3','PFI0','PFI1','PFI2','PFI3'},...
             'Whisker',4,'OutlierSize',1,'Symbol','r.');
else
    boxplot(a1,b1,...
             'Labels',{'K-S','R','D-R','PPA0','PPA2','PPA3','PPA4','PFI0','PFI1','PFI2','PFI3','PFI4'},...
             'Whisker',4,'OutlierSize',1,'Symbol','r.');
end
title('Panel A')
if case_nr==1
    axis([0.5 11.5 -5e-4 5e-4])
else
    axis([0.5 12.5 -2e-4 2e-4])
end
subplot(1,2,2)
if case_nr==1
    boxplot(a2,b2,...
             'Labels',{'K-S','R','D-R','PPA0','PPA1','PPA2','PPA3','PFI0','PFI1','PFI2','PFI3'},...
             'Whisker',4,'OutlierSize',1,'Symbol','r.');
else
    boxplot(a2,b2,...
             'Labels',{'K-S','R','D-R','PPA0','PPA2','PPA3','PPA4','PFI0','PFI1','PFI2','PFI3','PFI4'},...
             'Whisker',4,'OutlierSize',1,'Symbol','r.');
end
title('Panel B')   
if case_nr==1
    axis([0.5 11.5 -5e-4 5e-4])
else
    axis([0.5 12.5 -3e-4 3e-4])
end
clear a1 b1 a2 b2;


%--------------------------------------------------------------------------
% Boxplot of law of motion errors
%--------------------------------------------------------------------------
if case_nr==1
    a1=[Res.LoMErr_KS1;...
        Res.LoMErr_PPA0;Res.LoMErr_PPA1;Res.LoMErr_PPA2;Res.LoMErr_PPA3;...
        Res.LoMErr_PFI0;Res.LoMErr_PFI1;Res.LoMErr_PFI2;Res.LoMErr_PFI3];
    b1=[zeros(length(Res.LoMErr_KS1),1);ones(length(Res.LoMErr_PPA0),1);...
       2*ones(length(Res.LoMErr_PPA1),1);3*ones(length(Res.LoMErr_PPA2),1);4*ones(length(Res.LoMErr_PPA3),1);...
       5*ones(length(Res.LoMErr_PFI0),1);6*ones(length(Res.LoMErr_PFI1),1);...
       7*ones(length(Res.LoMErr_PFI2),1);8*ones(length(Res.LoMErr_PFI3),1)];
else
    a1=[Res.LoMErr_KS1;...
        Res.LoMErr_PPA0;Res.LoMErr_PPA2;Res.LoMErr_PPA3;Res.LoMErr_PPA4;...
        Res.LoMErr_PFI0;Res.LoMErr_PFI1;Res.LoMErr_PFI2;Res.LoMErr_PFI3;Res.LoMErr_PFI4];
    b1=[zeros(length(Res.LoMErr_KS1),1);ones(length(Res.LoMErr_PPA0),1);...
       2*ones(length(Res.LoMErr_PPA2),1);3*ones(length(Res.LoMErr_PPA3),1);...
       4*ones(length(Res.LoMErr_PPA4),1);5*ones(length(Res.LoMErr_PFI0),1);...
       6*ones(length(Res.LoMErr_PFI1),1);7*ones(length(Res.LoMErr_PFI2),1);...
       8*ones(length(Res.LoMErr_PFI3),1);9*ones(length(Res.LoMErr_PFI4),1)];
end
figure(4*(case_nr-1)+6)
hold on
plot([0 150],StaticParams.criter_k*ones(1,2),'k',[0 150],-StaticParams.criter_k*ones(1,2),'k')
if case_nr==1
    boxplot(a1,b1,...
             'Labels',{'K-S','PPA0','PPA1','PPA2','PPA3','PFI0','PFI1','PFI2','PFI3'},...
             'Whisker',4,'OutlierSize',1,'Symbol','r.');
else
    boxplot(a1,b1,...
         'Labels',{'K-S','PPA0','PPA2','PPA3','PPA4','PFI0','PFI1','PFI2','PFI3','PFI4'},...
         'Whisker',4,'OutlierSize',1,'Symbol','r.');
end
hold off
if case_nr==1
    axis([0.5 9.5 -5.1e-4 4.5e-4])
else
    axis([0.5 10.5 -1e-4 2e-4])
end


%--------------------------------------------------------------------------
% Stationary distributions 
%--------------------------------------------------------------------------
% preSol = load(strcat(folder,'/preSol_PPA.mat'));
Sol0 = load(strcat(folder,'/Sol0_PPA.mat'));
Sol1 = load(strcat(folder,'/Sol1_PPA.mat'));
Sol2 = load(strcat(folder,'/Sol2_PPA.mat'));
Sol3 = load(strcat(folder,'/Sol3_PPA.mat'));
if case_nr==2
    Sol4 = load(strcat(folder,'/Sol4_PPA.mat'));
end
KS = load(strcat(folder,'/KS1_Sol',num2str(case_nr),'.mat'));
KS2 = load(strcat(folder,'/KS2_Sol',num2str(case_nr),'.mat'));
KS3 = load(strcat(folder,'/KS3_Sol',num2str(case_nr),'.mat'));
if case_nr==1
    KS4 = load(strcat(folder,'/KS4_Sol',num2str(case_nr),'.mat'));
end
R = load(strcat(folder,'/R_Sol',num2str(case_nr),'.mat'));

p=[sum(StaticParams.p([1,3])),sum(StaticParams.p([2,4]))];

figure(4*(case_nr-1)+7)
hold on
if case_nr==2
    plot(StaticParams.kGrid,p*(Sol2.pdf-Sol0.pdf),'b')
    plot(StaticParams.kGrid,p*(Sol3.pdf-Sol2.pdf),'r')
    plot(StaticParams.kGrid,p*(Sol4.pdf-Sol3.pdf),'g')
else
    plot(StaticParams.kGrid,p*(Sol1.pdf-Sol0.pdf),'b')
    plot(StaticParams.kGrid,p*(Sol2.pdf-Sol1.pdf),'r')
    plot(StaticParams.kGrid,p*(Sol3.pdf-Sol2.pdf),'g')
end
plot([0 150],StaticParams.criter_k*ones(1,2),'k',[0 150],-StaticParams.criter_k*ones(1,2),'k')
hold off 
if case_nr==1
    axis([0 150 -1.2e-4 1.2e-4])
    legend('1^{st} order - 0^{th} order','2^{nd} order - 1^{st} order','3^{rd} order - 2^{nd} order', 'termination criterion','Location','northeast')
else
    axis([0 150 -1.2e-4 1.5e-4])
    legend('2^{nd} order - 0^{th} order', '3^{rd} order - 2^{nd} order', '4^{th} order - 3^{rd} order','termination criterion','Location','northeast')
end
xlabel('individual capital k');
ylabel('change in probability');

figure(4*(case_nr-1)+8)
subplot(1,2,1)
hold on
% plot(StaticParams.kGrid,(StaticParams.p'*[preSol.pdf_cond_b;preSol.pdf_cond_g]),'r','LineWidth',1);
plot(StaticParams.kGrid,p*KS.pdf,'k','LineWidth',1);
plot(StaticParams.kGrid,p*R.pdf,'g','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol0.pdf,'r','LineWidth',1);
if case_nr==1
    plot(StaticParams.kGrid,p*Sol1.pdf,'c','LineWidth',1.5);
    plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
    plot(StaticParams.kGrid,p*Sol3.pdf,'m--','LineWidth',1.5);
else
    plot(StaticParams.kGrid,p*Sol2.pdf,'c','LineWidth',1.5);
    plot(StaticParams.kGrid,p*Sol3.pdf,'b','LineWidth',1.5);
    plot(StaticParams.kGrid,p*Sol4.pdf,'m--','LineWidth',1.5);
end
hold off
if case_nr==1
    axis([0 100 0 0.0025])
    legend('K-S','R','PPA 0^{th} order','PPA 1^{st} order','PPA 2^{nd} order','PPA 3^{rd} order')
else
    axis([0 100 0 0.003])
    legend('K-S','R','PPA 0^{th} order','PPA 2^{nd} order','PPA 3^{rd} order','PPA 4^{th} order')
end
xlabel('individual capital k');
ylabel('probability');
title('Panel A')
subplot(1,2,2)
hold on
if case_nr==1
    plot(StaticParams.kGrid,p*Sol2.pdf,'k','LineWidth',1);
else
    plot(StaticParams.kGrid,p*Sol3.pdf,'k','LineWidth',1);
end
plot(StaticParams.kGrid,p*KS.pdf,'r','LineWidth',1);
plot(StaticParams.kGrid,p*KS2.pdf,'c','LineWidth',1.5);
plot(StaticParams.kGrid,p*KS3.pdf,'b','LineWidth',1.5);
if case_nr==1
    plot(StaticParams.kGrid,p*KS4.pdf,'m--','LineWidth',1.5);
end
hold off
if case_nr==1
    axis([0 100 0 0.0025])
    legend('PPA 2^{nd} order','K-S 1 moment','K-S 2 moments','K-S 3 moments','K-S 4 moments')
else
    axis([0 100 0 0.025])
    legend('PPA 3^{rd} order','K-S 1 moment','K-S 2 moments','K-S 3 moments')
end
xlabel('individual capital k');
ylabel('probability');
title('Panel B')
end

%--------------------------------------------------------------------------
% High vs. low unempl. benefit
%--------------------------------------------------------------------------

StaticParams = load(strcat('res_case1/StaticParams.mat'));
preSol1 = load(strcat('res_case1/preSol_PPA.mat'));
preSol2 = load(strcat('res_case2/preSol_PPA.mat'));
preSol3 = load(strcat('res_case3/preSol_PPA.mat'));
preSol4 = load(strcat('res_case4/preSol_PPA.mat'));
Sol1 = load(strcat('res_case1/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case2/Sol3_PPA.mat'));
Sol3 = load(strcat('res_case3/Sol2_PPA.mat'));
Sol4 = load(strcat('res_case4/Sol3_PPA.mat'));

figure(13)
hold on 
plot(StaticParams.kGrid,p*Sol3.pdf,'b','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol4.pdf,'g','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'m','LineWidth',1.5);
hold off
axis([0 100 0 0.003])
legend('\nu=5%','\nu=15%','\nu=40%','\nu=65%')
xlabel('individual capital k');
ylabel('probability');


figure(14)
subplot(1,2,1)
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol1.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol1.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol1.pdf,'g','LineWidth',1.5);
hold off
axis([0 130 0 0.003])
legend('no aggr. risk: bad state','no aggr. risk: good state','aggr. risk with both states')
xlabel('individual capital k');
ylabel('probability');
title('Panel A')%\nu=15%

subplot(1,2,2)
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol2.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol2.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol2.pdf,'g','LineWidth',1.5);
hold off
axis([0 140 0 0.0025])
xlabel('individual capital k');
ylabel('probability');
title('Panel B')%\nu=65%



%--------------------------------------------------------------------------
% Robustness Check
%--------------------------------------------------------------------------
% 1 benchmark model (see online appendix on model parameters)
% 2 unemployment benefit of 65%
% 3 unemployment benefit of 5%
% 4 unemployment benefit of 40%
% 5 subjective discount factor of 0.9
% 6 subjective discount factor of 0.999
% 7 risk aversion of 2
% 8 risk aversion of 4
% 9 production function elasticity of 0.2
% 10 production function elasticity of 0.45
% 11 depreciation rate of 1%
% 12 depreciation rate of 7.5%
% 13 unemployment benefit of 65%, subjective discount factor of 0.9
% 14 unemployment benefit of 65%, subjective discount factor of 0.999
% 15 unemployment benefit of 65%, risk aversion of 2
% 16 unemployment benefit of 65%, risk aversion of 4
% 17 unemployment benefit of 65%, production function elasticity of 0.2
% 18 unemployment benefit of 65%, production function elasticity of 0.45
% 19 unemployment benefit of 65%, depreciation rate of 1%
% 20 unemployment benefit of 65%, depreciation rate of 7.5%

figure(15)
subplot(4,3,1)
Sol1 = load(strcat('res_case5/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case13/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 20 0 0.021])
legend('\nu=15%','\nu=65%')
xlabel('individual capital k');
ylabel('probability');
title('\beta = 0.9')

subplot(4,3,2)
Sol1 = load(strcat('res_case1/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case2/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 100 0 0.003])
xlabel('individual capital k');
ylabel('probability');
title('Benchmark \beta = 0.99')

subplot(4,3,3)
Sol1 = load(strcat('res_case6/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case14/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 140 0 0.002])
xlabel('individual capital k');
ylabel('probability');
title('\beta = 0.999')

subplot(4,3,5)
Sol1 = load(strcat('res_case7/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case15/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 100 0 0.003])
xlabel('individual capital k');
ylabel('probability');
title('\gamma = 2')

subplot(4,3,4)
Sol1 = load(strcat('res_case1/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case2/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 100 0 0.003])
xlabel('individual capital k');
ylabel('probability');
title('Benchmark \gamma = 1')

subplot(4,3,6)
Sol1 = load(strcat('res_case8/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case16/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 100 0 0.003])
xlabel('individual capital k');
ylabel('probability');
title('\gamma = 4')

subplot(4,3,7)
Sol1 = load(strcat('res_case9/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case17/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 30 0 0.014])
xlabel('individual capital k');
ylabel('probability');
title('\alpha = 0.2')

subplot(4,3,8)
Sol1 = load(strcat('res_case1/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case2/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 100 0 0.003])
xlabel('individual capital k');
ylabel('probability');
title('Benchmark \alpha = 0.36')

subplot(4,3,9)
Sol1 = load(strcat('res_case10/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case18/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 300 0 0.001])
xlabel('individual capital k');
ylabel('probability');
title('\alpha = 0.45')

subplot(4,3,10)
Sol1 = load(strcat('res_case11/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case19/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 300 0 0.001])
xlabel('individual capital k');
ylabel('probability');
title('\delta = 0.01')

subplot(4,3,11)
Sol1 = load(strcat('res_case1/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case2/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 100 0 0.003])
xlabel('individual capital k');
ylabel('probability');
title('Benchmark \delta = 0.025')

subplot(4,3,12)
Sol1 = load(strcat('res_case12/Sol2_PPA.mat'));
Sol2 = load(strcat('res_case20/Sol3_PPA.mat'));
hold on 
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
hold off
axis([0 30 0 0.012])
xlabel('individual capital k');
ylabel('probability');
title('\delta = 0.075')

% 5 subjective discount factor of 0.9
% 6 subjective discount factor of 0.999
% 7 risk aversion of 2
% 8 risk aversion of 4
% 9 production function elasticity of 0.2
% 10 production function elasticity of 0.45
% 11 depreciation rate of 1%
% 12 depreciation rate of 7.5%
% 13 unemployment benefit of 65%, subjective discount factor of 0.9
% 14 unemployment benefit of 65%, subjective discount factor of 0.999
% 15 unemployment benefit of 65%, risk aversion of 2
% 16 unemployment benefit of 65%, risk aversion of 4
% 17 unemployment benefit of 65%, production function elasticity of 0.2
% 18 unemployment benefit of 65%, production function elasticity of 0.45
% 19 unemployment benefit of 65%, depreciation rate of 1%
% 20 unemployment benefit of 65%, depreciation rate of 7.5%

figure(15)
subplot(4,4,1)
Sol = load(strcat('res_case5/Sol3_PPA.mat'));
preSol = load(strcat('res_case5/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 20 0 0.021])
legend('no aggr. risk: bad state','no aggr. risk: good state','aggr. risk with both states')
xlabel('individual capital k');
ylabel('probability');
title('\beta = 0.9, \nu=15%')

subplot(4,4,2)
Sol = load(strcat('res_case13/Sol3_PPA.mat'));
preSol = load(strcat('res_case13/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 20 0 0.021])
xlabel('individual capital k');
ylabel('probability');
title('\beta = 0.9, \nu=65%')

subplot(4,4,3)
Sol = load(strcat('res_case6/Sol3_PPA.mat'));
preSol = load(strcat('res_case6/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 160 0 0.0031])
xlabel('individual capital k');
ylabel('probability');
title('\beta = 0.999, \nu=15%')

subplot(4,4,4)
Sol = load(strcat('res_case14/Sol3_PPA.mat'));
preSol = load(strcat('res_case14/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 160 0 0.0031])
xlabel('individual capital k');
ylabel('probability');
title('\beta = 0.999, \nu=65%')

subplot(4,4,5)
Sol = load(strcat('res_case7/Sol3_PPA.mat'));
preSol = load(strcat('res_case7/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 140 0 0.0031])
xlabel('individual capital k');
ylabel('probability');
title('\gamma = 2, \nu=15%')

subplot(4,4,6)
Sol = load(strcat('res_case15/Sol3_PPA.mat'));
preSol = load(strcat('res_case15/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 140 0 0.0031])
xlabel('individual capital k');
ylabel('probability');
title('\gamma = 2, \nu=65%')

subplot(4,4,7)
Sol = load(strcat('res_case8/Sol3_PPA.mat'));
preSol = load(strcat('res_case8/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 140 0 0.0035])
xlabel('individual capital k');
ylabel('probability');
title('\gamma = 4, \nu=15%')

subplot(4,4,8)
Sol = load(strcat('res_case16/Sol3_PPA.mat'));
preSol = load(strcat('res_case16/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 140 0 0.0035])
xlabel('individual capital k');
ylabel('probability');
title('\gamma = 4, \nu=65%')

subplot(4,4,9)
Sol = load(strcat('res_case9/Sol3_PPA.mat'));
preSol = load(strcat('res_case9/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 30 0 0.013])
xlabel('individual capital k');
ylabel('probability');
title('\alpha = 0.2, \nu=15%')

subplot(4,4,10)
Sol = load(strcat('res_case17/Sol3_PPA.mat'));
preSol = load(strcat('res_case17/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 30 0 0.013])
xlabel('individual capital k');
ylabel('probability');
title('\alpha = 0.2, \nu=65%')

subplot(4,4,11)
Sol = load(strcat('res_case10/Sol3_PPA.mat'));
preSol = load(strcat('res_case10/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 300 0 0.001])
xlabel('individual capital k');
ylabel('probability');
title('\alpha = 0.45, \nu=15%')

subplot(4,4,12)
Sol = load(strcat('res_case18/Sol3_PPA.mat'));
preSol = load(strcat('res_case18/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 300 0 0.001])
xlabel('individual capital k');
ylabel('probability');
title('\alpha = 0.45, \nu=65%')

subplot(4,4,13)
Sol = load(strcat('res_case11/Sol3_PPA.mat'));
preSol = load(strcat('res_case11/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 300 0 0.001])
xlabel('individual capital k');
ylabel('probability');
title('\delta = 1%, \nu=15%')

subplot(4,4,14)
Sol = load(strcat('res_case19/Sol3_PPA.mat'));
preSol = load(strcat('res_case19/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 300 0 0.001])
xlabel('individual capital k');
ylabel('probability');
title('\delta = 1%, \nu=65%')

subplot(4,4,15)
Sol = load(strcat('res_case12/Sol3_PPA.mat'));
preSol = load(strcat('res_case12/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 30 0 0.019])
xlabel('individual capital k');
ylabel('probability');
title('\delta = 7.5%, \nu=15%')

subplot(4,4,16)
Sol = load(strcat('res_case20/Sol3_PPA.mat'));
preSol = load(strcat('res_case20/preSol_PPA.mat'));
hold on 
plot(StaticParams.kGrid,2*StaticParams.p(1:2)'*preSol.pdf_cond_b,'r')
plot(StaticParams.kGrid,2*StaticParams.p(3:4)'*preSol.pdf_cond_g,'b')
plot(StaticParams.kGrid,p*Sol.pdf,'g','LineWidth',1.5);
hold off
axis([0 30 0 0.019])
xlabel('individual capital k');
ylabel('probability');
title('\delta = 7.5%, \nu=65%')