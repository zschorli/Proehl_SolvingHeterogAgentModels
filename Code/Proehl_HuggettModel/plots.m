% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This file produces the plots in the paper.
%__________________________________________________________________________
clc;
clear all;

case_nr = 1;
folder = strcat('res_Huggett_case',num2str(case_nr));
StaticParams = load(strcat(folder,'/StaticParams.mat'));

%--------------------------------------------------------------------------
% Base distribution
%--------------------------------------------------------------------------
Poly =load(strcat(folder,'/Poly3_PPA.mat'));
figure(1)
subplot(1,2,1)
plot(sort([StaticParams.kGrid,StaticParams.kGrid(1:end-1)]),Poly.pdf{2}(1,sort([1:length(Poly.pdf{2}),2:length(Poly.pdf{2})])),'r');
axis([StaticParams.k_min StaticParams.k_max 0 0.62])
title('\xi^a_1');
xlabel('individual asset allocation a');
ylabel('probability');

subplot(1,2,2)
plot(sort([StaticParams.kGrid,StaticParams.kGrid(1:end-1)]),Poly.pdf{3}(1,sort([1:length(Poly.pdf{2}),2:length(Poly.pdf{2})])),'b');
axis([StaticParams.k_min StaticParams.k_max 0 0.08])
title('\xi^a_2');
xlabel('individual asset allocation a');
ylabel('probability');


%--------------------------------------------------------------------------
% Orthogonal Polynomials
%--------------------------------------------------------------------------
figure(2)
subplot(1,3,1)
hold on
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{2}(1,:)*Poly.xi{2},'b')
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{2}(2,:)*Poly.xi{2},'g')
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{2}(3,:)*Poly.xi{2},'r')
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{2}(4,:)*Poly.xi{2},'c')
axis([StaticParams.k_min 4 -2.5 2.2])
legend('\Phi^{k_1}_0','\Phi^{k_1}_1','\Phi^{k_1}_2','\Phi^{k_1}_3','Location','SouthEast');
xlabel('basic random variable \xi^k_1');
ylabel('function value of the polynomial');
title('Panel A')
subplot(1,3,2)
hold on
plot(Poly.xi{3}(2,:),Poly.phi_coeffs{3}(1,:)*Poly.xi{3},'b')
plot(Poly.xi{3}(2,:),Poly.phi_coeffs{3}(2,:)*Poly.xi{3},'g')
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{3}(3,:)*Poly.xi{3},'r')
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{3}(4,:)*Poly.xi{3},'c')
axis([StaticParams.k_min 4 -2.5 2.2])
legend('\Phi^{k_2}_0','\Phi^{k_2}_1','\Phi^{k_2}_2','\Phi^{k_2}_3','Location','SouthEast');
xlabel('basic random variable \xi^k_2');
ylabel('function value of the polynomial');
title('Panel B')
subplot(1,3,3)
hold on
plot(Poly.xi{1}(2,[1,3]),Poly.phi_coeffs{1}(1,:)*Poly.xi{1}(:,[1,3]),'b*')
plot(Poly.xi{1}(2,[1,3]),Poly.phi_coeffs{1}(2,:)*Poly.xi{1}(:,[1,3]),'g*')
plot(Poly.xi{1}(2,[1,3]),Poly.phi_coeffs{1}(3,:)*Poly.xi{1}(:,[1,3]),'r*')
hold off
axis([0.5 3.5 -1.2 1.2])
legend('\Phi^z_0','\Phi^z_1','\Phi^z_2','Location','SouthEast');
xlabel('basic random variable \xi^z');
ylabel('function value of the polynomial');
title('Panel C')


for case_nr = 1:2
switch case_nr
    case 1
        folder = 'res_Huggett_case1';
    case 2
        folder = 'res_Huggett_case2';
end
StaticParams = load(strcat(folder,'/StaticParams.mat'));
%--------------------------------------------------------------------------
% Boxplot of Euler equation errors
%--------------------------------------------------------------------------
Res = load(strcat(folder,'/Res.mat'));

a1=[Res.stErr_PPA1;Res.stErr_PPA2;Res.stErr_PPA3;Res.stErr_PPA4;...
    Res.stErr_PFI1;Res.stErr_PFI2;Res.stErr_PFI3;Res.stErr_PFI4];
b1=[zeros(length(Res.stErr_PPA1),1);ones(length(Res.stErr_PPA2),1);...
    2*ones(length(Res.stErr_PPA3),1);3*ones(length(Res.stErr_PPA4),1);...
    4*ones(length(Res.stErr_PFI1),1);5*ones(length(Res.stErr_PFI2),1);...
    6*ones(length(Res.stErr_PFI3),1);7*ones(length(Res.stErr_PFI4),1)];
a2=[Res.dynErr_PPA1;Res.dynErr_PPA2;Res.dynErr_PPA3;Res.dynErr_PPA4;...
    Res.dynErr_PFI1;Res.dynErr_PFI2;Res.dynErr_PFI3;Res.dynErr_PFI4];
b2=[zeros(length(Res.dynErr_PPA1),1);ones(length(Res.dynErr_PPA2),1);...
    2*ones(length(Res.dynErr_PPA3),1);3*ones(length(Res.dynErr_PPA4),1);...
    4*ones(length(Res.dynErr_PFI1),1);5*ones(length(Res.dynErr_PFI2),1);...
    6*ones(length(Res.dynErr_PFI3),1);7*ones(length(Res.dynErr_PFI4),1)];  

figure(7*(case_nr-1)+3)
subplot(1,2,1)
boxplot(abs(a1),b1,...
         'Labels',{'PPA1','PPA2','PPA3','PPA4','PFI1','PFI2','PFI3','PFI4'},...
         'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
if case_nr==1
    axis([0.5 8.5 -0.005 0.15])
else
    axis([0.5 8.5 -0.005 0.06])
end
subplot(1,2,2)
boxplot(abs(a2),b2,...
         'Labels',{'PPA1','PPA2','PPA3','PPA4','PFI1','PFI2','PFI3','PFI4'},...
         'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel B')   
if case_nr==1
    axis([0.5 8.5 -0.005 0.3])
else
    axis([0.5 8.5 -0.005 0.07])
end

figure(7*(case_nr-1)+4)
subplot(1,2,1)
boxplot(a1,b1,...
        'Labels',{'PPA1','PPA2','PPA3','PPA4','PFI1','PFI2','PFI3','PFI4'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
if case_nr==1
    axis([0.5 8.5 -0.06 0.18])
else
    axis([0.5 8.5 -0.06 0.04])
end
subplot(1,2,2)
boxplot(a2,b2,...
        'Labels',{'PPA1','PPA2','PPA3','PPA4','PFI1','PFI2','PFI3','PFI4'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel B')   
if case_nr==1
    axis([0.5 8.5 -0.15 0.37])
else
    axis([0.5 8.5 -0.081 0.08])
end
clear a1 b1 a2 b2;


%--------------------------------------------------------------------------
% Boxplot of price errors
%--------------------------------------------------------------------------
a1=[Res.stErr_PPA1_pr;Res.stErr_PPA2_pr;Res.stErr_PPA3_pr;Res.stErr_PPA4_pr;...
    Res.stErr_PFI1_pr;Res.stErr_PFI2_pr;Res.stErr_PFI3_pr;Res.stErr_PFI4_pr];
b1=[zeros(length(Res.stErr_PPA1_pr),1);ones(length(Res.stErr_PPA2_pr),1);...
    2*ones(length(Res.stErr_PPA3_pr),1);3*ones(length(Res.stErr_PPA4_pr),1);...
    4*ones(length(Res.stErr_PFI1_pr),1);5*ones(length(Res.stErr_PFI2_pr),1);...
    6*ones(length(Res.stErr_PFI3_pr),1);7*ones(length(Res.stErr_PFI4_pr),1)];
a2=[Res.dynErr_PPA1_pr;Res.dynErr_PPA2_pr;Res.dynErr_PPA3_pr;Res.dynErr_PPA4_pr;...
    Res.dynErr_PFI1_pr;Res.dynErr_PFI2_pr;Res.dynErr_PFI3_pr;Res.dynErr_PFI4_pr];
b2=[zeros(length(Res.dynErr_PPA1_pr),1);ones(length(Res.dynErr_PPA2_pr),1);...
    2*ones(length(Res.dynErr_PPA3_pr),1);3*ones(length(Res.dynErr_PPA4_pr),1);...
    4*ones(length(Res.dynErr_PFI1_pr),1);5*ones(length(Res.dynErr_PFI2_pr),1);...
    6*ones(length(Res.dynErr_PFI3_pr),1);7*ones(length(Res.dynErr_PFI4_pr),1)];  

figure(7*(case_nr-1)+5)
subplot(1,2,1)
boxplot(abs(a1),b1,...
         'Labels',{'PPA1','PPA2','PPA3','PPA4','PFI1','PFI2','PFI3','PFI4'},...
         'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
if case_nr==1
    axis([0.5 8.5 -1e-2 0.3])
else
    axis([0.5 8.5 -1e-2 0.12])
end
subplot(1,2,2)
boxplot(abs(a2),b2,...
         'Labels',{'PPA1','PPA2','PPA3','PPA4','PFI1','PFI2','PFI3','PFI4'},...
         'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel B')   
if case_nr==1
    axis([0.5 8.5 -1e-2 0.4])
else
    axis([0.5 8.5 -1e-2 0.18])
end


figure(7*(case_nr-1)+6)
subplot(1,2,1)
boxplot(a1,b1,...
        'Labels',{'PPA1','PPA2','PPA3','PPA4','PFI1','PFI2','PFI3','PFI4'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
if case_nr==1
    axis([0.5 8.5 -0.02 0.3])
else
    axis([0.5 8.5 -0.07 0.12])
end
subplot(1,2,2)
boxplot(a2,b2,...
        'Labels',{'PPA1','PPA2','PPA3','PPA4','PFI1','PFI2','PFI3','PFI4'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel B')   
if case_nr==1
    axis([0.5 8.5 0 0.4])
else
    axis([0.5 8.5 -0.1 0.19])
end


%--------------------------------------------------------------------------
% Boxplot of Law of motion errors
%--------------------------------------------------------------------------
a1=[squeeze(Res.K(1,1,:));squeeze(Res.K(1,2,:));squeeze(Res.K(1,3,:));squeeze(Res.K(1,4,:));...
    squeeze(Res.K(2,1,:));squeeze(Res.K(2,2,:));squeeze(Res.K(2,3,:));squeeze(Res.K(2,4,:))];
b1=[zeros(length(squeeze(Res.K(1,1,:))),1);ones(length(squeeze(Res.K(1,1,:))),1);...
   2*ones(length(squeeze(Res.K(1,1,:))),1);3*ones(length(squeeze(Res.K(1,1,:))),1);...
   4*ones(length(squeeze(Res.K(1,1,:))),1);5*ones(length(squeeze(Res.K(1,1,:))),1);...
   6*ones(length(squeeze(Res.K(1,1,:))),1);7*ones(length(squeeze(Res.K(1,1,:))),1)];

figure(7*(case_nr-1)+7)
hold on
boxplot(a1,b1,...
        'Labels',{'PPA1','PPA2','PPA3','PPA4','PFI1','PFI2','PFI3','PFI4'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
hold off
if case_nr==1
    axis([0.5 8.5 -0.1 0.045])
else
    axis([0.5 8.5 -0.02 0.05])
end
clear a1 b1 a2 b2;


%--------------------------------------------------------------------------
% Stationary distributions 
%--------------------------------------------------------------------------
% preSol = load(strcat(folder,'/preSol_PPA.mat'));
Sol1 = load(strcat(folder,'/Sol1_PPA.mat'));
Sol2 = load(strcat(folder,'/Sol2_PPA.mat'));
Sol3 = load(strcat(folder,'/Sol3_PPA.mat'));
Sol4 = load(strcat(folder,'/Sol4_PPA.mat'));

p=[sum(StaticParams.p([1,3])),sum(StaticParams.p([2,4]))];

figure(7*(case_nr-1)+8)
hold on
plot(StaticParams.kGrid,p*(Sol2.pdf-Sol1.pdf),'b')
plot(StaticParams.kGrid,p*(Sol3.pdf-Sol2.pdf),'r')
plot(StaticParams.kGrid,p*(Sol4.pdf-Sol3.pdf),'g')
plot([StaticParams.k_min StaticParams.k_max],StaticParams.criter_k*ones(1,2),'k',[StaticParams.k_min StaticParams.k_max],-StaticParams.criter_k*ones(1,2),'k')
hold off 
axis([StaticParams.k_min StaticParams.k_max -1.8e-4 1.8e-4])
legend('2^{nd} order - 1^{st} order','3^{rd} order - 2^{nd} order','4^{th} order - 3^{rd} order','termination criterion','Location','northeast')%'4^{th} order - 3^{rd} order', 
xlabel('individual asset allocation a');
ylabel('change in probability');

figure(7*(case_nr-1)+9)
subplot(1,2,1)
hold on
plot(StaticParams.kGrid,p*Sol1.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'c','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol3.pdf,'b','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol4.pdf,'m--','LineWidth',1.5);
hold off
if case_nr==1
    axis([StaticParams.k_min 3 0 0.045])
else
    axis([StaticParams.k_min 1 0 0.12])
end
legend('PPA 1^{st} order','PPA 2^{nd} order','PPA 3^{rd} order','PPA 4^{th} order')%
xlabel('individual asset allocation a');
ylabel('probability');
title('Panel A')
subplot(1,2,2)
hold on
x = sort([StaticParams.kGrid(1:end),StaticParams.kGrid(1:end-1)+1e-4]);
y = p*Sol1.pdf;
plot(x,y(sort([2:end,1:end])),'r','LineWidth',1.5);
y = p*Sol2.pdf;
plot(x,y(sort([2:end,1:end])),'c','LineWidth',1.5);
y = p*Sol3.pdf;
plot(x,y(sort([2:end,1:end])),'b','LineWidth',1.5);
y = p*Sol4.pdf;
plot(x,y(sort([2:end,1:end])),'m--','LineWidth',1.5);
hold off
if case_nr==1
    axis([StaticParams.k_min -1 0 0.045])
else
    axis([StaticParams.k_min -0.3 0 0.12])
end
xlabel('individual asset allocation a');
ylabel('probability');
title('Panel B')
end

%--------------------------------------------------------------------------
% Price Distribution
%--------------------------------------------------------------------------
Res1 = load('res_Huggett_case1/Res.mat');
[a1,b1]=hist(squeeze(Res1.pr(1,4,:)),100);
a1=a1/sum(a1);
Res2 = load('res_Huggett_case2/Res.mat');
[a2,b2]=hist(squeeze(Res2.pr(1,4,:)),100);
a2=a2/sum(a2);
Res3 = load('res_Huggett_case3/Res.mat');
[a3,b3]=hist(squeeze(Res3.pr(1,4,:)),100);
a3=a3/sum(a3);
Res4 = load('res_Huggett_case4/Res.mat');
[a4,b4]=hist(squeeze(Res4.pr(1,4,:)),100);
a4=a4/sum(a4);
Res5 = load('res_Huggett_case5/Res.mat');
[a5,b5]=hist(squeeze(Res5.pr(1,3,:)),100);
a5=a5/sum(a5);


figure(17)
hold on
h4=plot([b4(1),sort([b4,b4(2:end)-1e-4]),b4(end)],[0,a4(sort([1:end,1:end-1])),0],'g','LineWidth',1.5);
h2=plot([b3(1),sort([b3,b3(2:end)-1e-4]),b3(end)],[0,a3(sort([1:end,1:end-1])),0],'c','LineWidth',1.5);
h3=plot([b1(1),sort([b1,b1(2:end)-1e-4]),b1(end)],[0,a1(sort([1:end,1:end-1])),0],'r');
h1=plot([b2(1),sort([b2,b2(2:end)-1e-4]),b2(end)],[0,a2(sort([1:end,1:end-1])),0],'b');
h5=plot([b5(1),sort([b5,b5(2:end)-1e-4]),b5(end)],[0,a5(sort([1:end,1:end-1])),0],'m');
hold off
xlabel('prices');
ylabel('probability');
legend([h1,h2,h3,h4,h5],{'$\bar{a}=-0.5$','$\bar{a}=-1$','$\bar{a}=-1.5$','$\bar{a}=-2$','$\bar{a}=-2.5$'},'Interpreter','latex','Location','Northwest')


%--------------------------------------------------------------------------
% Robustness of the Price Distribution
%--------------------------------------------------------------------------
figure(18)
Res6 = load('res_Huggett_case6/Res.mat');
[a1,b1]=hist(squeeze(Res6.pr(1,3,:)),100);
a1=a1/sum(a1);
Res7 = load('res_Huggett_case10/Res.mat');
[a2,b2]=hist(squeeze(Res7.pr(1,4,:)),100);
a2=a2/sum(a2);
subplot(2,3,1)
hold on
plot([b2(1),sort([b2,b2(2:end)-1e-4]),b2(end)],[0,a2(sort([1:end,1:end-1])),0],'b');
plot([b1(1),sort([b1,b1(2:end)-1e-4]),b1(end)],[0,a1(sort([1:end,1:end-1])),0],'r');
hold off
xlabel('prices');
ylabel('probability');
legend({'$\bar{a}=0.5$','$\bar{a}=1.5$'},'Interpreter','latex','Location','Northwest')
title('\beta=0.9')

Res6 = load('res_Huggett_case1/Res.mat');
[a1,b1]=hist(squeeze(Res6.pr(1,4,:)),100);
a1=a1/sum(a1);
Res7 = load('res_Huggett_case2/Res.mat');
[a2,b2]=hist(squeeze(Res7.pr(1,4,:)),100);
a2=a2/sum(a2);
subplot(2,3,2)
hold on
plot([b2(1),sort([b2,b2(2:end)-1e-4]),b2(end)],[0,a2(sort([1:end,1:end-1])),0],'b');
plot([b1(1),sort([b1,b1(2:end)-1e-4]),b1(end)],[0,a1(sort([1:end,1:end-1])),0],'r');
hold off
xlabel('prices');
ylabel('probability');
title('Benchmark \beta=0.99')

Res6 = load('res_Huggett_case7/Res.mat');
[a1,b1]=hist(squeeze(Res6.pr(1,4,:)),100);
a1=a1/sum(a1);
Res7 = load('res_Huggett_case11/Res.mat');
[a2,b2]=hist(squeeze(Res7.pr(1,4,:)),100);
a2=a2/sum(a2);
subplot(2,3,3)
hold on
plot([b2(1),sort([b2,b2(2:end)-1e-4]),b2(end)],[0,a2(sort([1:end,1:end-1])),0],'b');
plot([b1(1),sort([b1,b1(2:end)-1e-4]),b1(end)],[0,a1(sort([1:end,1:end-1])),0],'r');
hold off
xlabel('prices');
ylabel('probability');
title('\beta=0.999')

Res6 = load('res_Huggett_case8/Res.mat');
[a1,b1]=hist(squeeze(Res6.pr(1,4,:)),100);
a1=a1/sum(a1);
Res7 = load('res_Huggett_case12/Res.mat');
[a2,b2]=hist(squeeze(Res7.pr(1,4,:)),100);
a2=a2/sum(a2);
subplot(2,3,4)
hold on
plot([b2(1),sort([b2,b2(2:end)-1e-4]),b2(end)],[0,a2(sort([1:end,1:end-1])),0],'b');
plot([b1(1),sort([b1,b1(2:end)-1e-4]),b1(end)],[0,a1(sort([1:end,1:end-1])),0],'r');
hold off
xlabel('prices');
ylabel('probability');
title('\gamma=1')

Res6 = load('res_Huggett_case1/Res.mat');
[a1,b1]=hist(squeeze(Res6.pr(1,4,:)),100);
a1=a1/sum(a1);
Res7 = load('res_Huggett_case2/Res.mat');
[a2,b2]=hist(squeeze(Res7.pr(1,4,:)),100);
a2=a2/sum(a2);
subplot(2,3,5)
hold on
plot([b2(1),sort([b2,b2(2:end)-1e-4]),b2(end)],[0,a2(sort([1:end,1:end-1])),0],'b');
plot([b1(1),sort([b1,b1(2:end)-1e-4]),b1(end)],[0,a1(sort([1:end,1:end-1])),0],'r');
hold off
xlabel('prices');
ylabel('probability');
title('Benchmark \gamma=2')

Res6 = load('res_Huggett_case9/Res.mat');
[a1,b1]=hist(squeeze(Res6.pr(1,3,:)),100);
a1=a1/sum(a1);
Res7 = load('res_Huggett_case13/Res.mat');
[a2,b2]=hist(squeeze(Res7.pr(1,4,:)),100);
a2=a2/sum(a2);
subplot(2,3,6)
hold on
plot([b2(1),sort([b2,b2(2:end)-1e-4]),b2(end)],[0,a2(sort([1:end,1:end-1])),0],'b');
plot([b1(1),sort([b1,b1(2:end)-1e-4]),b1(end)],[0,a1(sort([1:end,1:end-1])),0],'r');
hold off
xlabel('prices');
ylabel('probability');
title('\gamma=3')
