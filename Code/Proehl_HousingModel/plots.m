% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
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
Poly =load(strcat(folder,'/Poly3_PFI.mat'));
figure(1)
subplot(1,2,1)
plot(sort([StaticParams.kGrid,StaticParams.kGrid(1:end-1)]),Poly.pdf{2}(1,sort([1:length(Poly.pdf{2}),2:length(Poly.pdf{2})])),'r');
axis([StaticParams.k_min 10 0 0.62])
title('\xi^a_1');
xlabel('individual net worth a');
ylabel('probability');

subplot(1,2,2)
plot(sort([StaticParams.kGrid,StaticParams.kGrid(1:end-1)]),Poly.pdf{3}(1,sort([1:length(Poly.pdf{2}),2:length(Poly.pdf{2})])),'b');
axis([StaticParams.k_min 10 0 0.02])
title('\xi^a_2');
xlabel('individual net worth a');
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
axis([StaticParams.k_min 6 -7 7])
legend('\Phi^{a_1}_0','\Phi^{a_1}_1','\Phi^{a_1}_2','\Phi^{a_1}_3','Location','SouthEast');
xlabel('basic random variable \xi^a_1');
ylabel('function value of the polynomial');
title('Panel A')
subplot(1,3,2)
hold on
plot(Poly.xi{3}(2,:),Poly.phi_coeffs{3}(1,:)*Poly.xi{3},'b')
plot(Poly.xi{3}(2,:),Poly.phi_coeffs{3}(2,:)*Poly.xi{3},'g')
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{3}(3,:)*Poly.xi{3},'r')
plot(Poly.xi{2}(2,:),Poly.phi_coeffs{3}(4,:)*Poly.xi{3},'c')
axis([StaticParams.k_min 6 -7 7])
legend('\Phi^{a_2}_0','\Phi^{a_2}_1','\Phi^{a_2}_2','\Phi^{a_2}_3','Location','SouthEast');
xlabel('basic random variable \xi^a_2');
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


folder = 'res_Huggett_case1';
StaticParams = load(strcat(folder,'/StaticParams.mat'));
%--------------------------------------------------------------------------
% Boxplot of Euler equation errors
%--------------------------------------------------------------------------
Res = load(strcat(folder,'/Res.mat'));

a1=[Res.stErr_PFI1;Res.stErr_PFI2;Res.stErr_PFI3;Res.stErr_PFI4];
b1=[zeros(length(Res.stErr_PFI1),1);ones(length(Res.stErr_PFI2),1);...
    2*ones(length(Res.stErr_PFI3),1);3*ones(length(Res.stErr_PFI4),1)];
a2=[Res.dynErr_PFI1;Res.dynErr_PFI2;Res.dynErr_PFI3;Res.dynErr_PFI4];
b2=[zeros(length(Res.dynErr_PFI1),1);ones(length(Res.dynErr_PFI2),1);...
    2*ones(length(Res.dynErr_PFI3),1);3*ones(length(Res.dynErr_PFI4),1)];  

figure(3)
subplot(1,2,1)
boxplot(a1,b1,...
        'Labels',{'PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
% axis([0.5 8.5 -0.06 0.18])
subplot(1,2,2)
boxplot(a2,b2,...
        'Labels',{'PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel B')   
% axis([0.5 8.5 -0.15 0.37])
clear a1 b1 a2 b2;


%--------------------------------------------------------------------------
% Boxplot of price errors
%--------------------------------------------------------------------------
a1=[Res.stErr_PFI1_q;Res.stErr_PFI2_q;Res.stErr_PFI3_q;Res.stErr_PFI4_q];
b1=[zeros(length(Res.stErr_PFI1_q),1);ones(length(Res.stErr_PFI2_q),1);...
    2*ones(length(Res.stErr_PFI3_q),1);3*ones(length(Res.stErr_PFI4_q),1)];
a2=[Res.dynErr_PFI1_q;Res.dynErr_PFI2_q;Res.dynErr_PFI3_q;Res.dynErr_PFI4_q];
b2=[zeros(length(Res.dynErr_PFI1_q),1);ones(length(Res.dynErr_PFI2_q),1);...
    2*ones(length(Res.dynErr_PFI3_q),1);3*ones(length(Res.dynErr_PFI4_q),1)];  

figure(4)
subplot(1,2,1)
boxplot(a1,b1,...
        'Labels',{'PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
% axis([0.5 8.5 -0.02 0.3])
subplot(1,2,2)
boxplot(a2,b2,...
        'Labels',{'PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel B')   
% axis([0.5 8.5 0 0.4])

a1=[Res.stErr_PFI1_H;Res.stErr_PFI2_H;Res.stErr_PFI3_H;Res.stErr_PFI4_H];
b1=[zeros(length(Res.stErr_PFI1_H),1);ones(length(Res.stErr_PFI2_H),1);...
    2*ones(length(Res.stErr_PFI3_H),1);3*ones(length(Res.stErr_PFI4_H),1)];
a2=[Res.dynErr_PFI1_H;Res.dynErr_PFI2_H;Res.dynErr_PFI3_H;Res.dynErr_PFI4_H];
b2=[zeros(length(Res.dynErr_PFI1_H),1);ones(length(Res.dynErr_PFI2_H),1);...
    2*ones(length(Res.dynErr_PFI3_H),1);3*ones(length(Res.dynErr_PFI4_H),1)];  

figure(5)
subplot(1,2,1)
boxplot(a1,b1,...
        'Labels',{'PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
% axis([0.5 8.5 -0.02 0.3])
subplot(1,2,2)
boxplot(a2,b2,...
        'Labels',{'PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel B')   
% axis([0.5 8.5 0 0.4])



%--------------------------------------------------------------------------
% Stationary distributions 
%--------------------------------------------------------------------------
Sol0 = load(strcat(folder,'/Sol0_PFI.mat'));
Sol1 = load(strcat(folder,'/Sol1_PFI.mat'));
Sol2 = load(strcat(folder,'/Sol2_PFI.mat'));
Sol3 = load(strcat(folder,'/Sol3_PFI.mat'));

p=[sum(StaticParams.p([1,3])),sum(StaticParams.p([2,4]))];

figure(6)
subplot(1,2,1)
hold on
plot(StaticParams.kGrid,p*(Sol1.pdf-Sol0.pdf),'b')
plot(StaticParams.kGrid,p*(Sol2.pdf-Sol1.pdf),'r')
plot(StaticParams.kGrid,p*(Sol3.pdf-Sol2.pdf),'g')
plot([StaticParams.k_min StaticParams.k_max],StaticParams.criter_k*ones(1,2),'k',[StaticParams.k_min StaticParams.k_max],-StaticParams.criter_k*ones(1,2),'k')
hold off 
axis([StaticParams.k_min 30 -6e-4 6e-4])
legend('1^{st} order - 0^{th} order','2^{nd} order - 1^{st} order','3^{rd} order - 2^{nd} order','termination criterion','Location','southeast')
xlabel('individual net worth a');
ylabel('change in probability');
subplot(1,2,2)
hold on
plot(StaticParams.kGrid,p*(Sol1.pdf-Sol0.pdf),'b')
plot(StaticParams.kGrid,p*(Sol2.pdf-Sol1.pdf),'r')
plot(StaticParams.kGrid,p*(Sol3.pdf-Sol2.pdf),'g')
plot([StaticParams.k_min StaticParams.k_max],StaticParams.criter_k*ones(1,2),'k',[StaticParams.k_min StaticParams.k_max],-StaticParams.criter_k*ones(1,2),'k')
hold off 
axis([StaticParams.k_min 1.5 -1.8e-3 7e-4])
xlabel('individual net worth a');
ylabel('change in probability');


figure(7)
subplot(1,2,1)
hold on
plot(StaticParams.kGrid,p*Sol0.pdf,'r','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol1.pdf,'c','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol2.pdf,'b','LineWidth',1.5);
plot(StaticParams.kGrid,p*Sol3.pdf,'m--','LineWidth',1.5);
hold off
axis([StaticParams.k_min 15 0 0.009])
legend('PFI 0^{th} order','PFI 1^{st} order','PFI 2^{nd} order','PFI 3^{rd} order')
xlabel('individual net worth a');
ylabel('probability');
title('Panel A')
subplot(1,2,2)
hold on
x = sort([StaticParams.kGrid(1:end),StaticParams.kGrid(1:end-1)+1e-4]);
y = p*Sol0.pdf;
plot(x,y(sort([2:end,1:end])),'r','LineWidth',1.5);
y = p*Sol1.pdf;
plot(x,y(sort([2:end,1:end])),'c','LineWidth',1.5);
y = p*Sol2.pdf;
plot(x,y(sort([2:end,1:end])),'b','LineWidth',1.5);
y = p*Sol3.pdf;
plot(x,y(sort([2:end,1:end])),'m--','LineWidth',1.5);
hold off
axis([StaticParams.k_min 0.08 0.2 0.26])
xlabel('individual net worth a');
ylabel('probability');
title('Panel B')


%--------------------------------------------------------------------------
% Price Distribution
%--------------------------------------------------------------------------
[a0,b0]=hist(squeeze(Res.q(end,Res.z_ag==0)),75);
[a1,b1]=hist(squeeze(Res.q(end,Res.z_ag==1)),75);

[c0,d0]=hist(squeeze(Res.r(end,Res.z_ag==0)),75);
[c1,d1]=hist(squeeze(Res.r(end,Res.z_ag==1)),75);

figure(8)
subplot(1,2,1)
hold on
h0=plot([b0(1),sort([b0,b0(2:end)-1e-4]),b0(end)],[0,a0(sort([1:end,1:end-1])),0],'r');
h1=plot([b1(1),sort([b1,b1(2:end)-1e-4]),b1(end)],[0,a1(sort([1:end,1:end-1])),0],'g');
hold off
xlabel('house price');
ylabel('probability');
legend([h0,h1],{'recessions','booms'},'Interpreter','latex','Location','Northwest')
subplot(1,2,2)
hold on
h0=plot([d0(1),sort([d0,d0(2:end)-1e-4]),d0(end)],[0,c0(sort([1:end,1:end-1])),0],'r');
h1=plot([d1(1),sort([d1,d1(2:end)-1e-4]),d1(end)],[0,c1(sort([1:end,1:end-1])),0],'g');
hold off
xlabel('interest rate');
ylabel('probability');
