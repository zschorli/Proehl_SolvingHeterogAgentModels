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

folder = strcat('res_Housing');
StaticParams = load(strcat(folder,'/StaticParams.mat'));

%--------------------------------------------------------------------------
% Base distribution
%--------------------------------------------------------------------------
Poly =load(strcat(folder,'/Poly3_PFI.mat'));
figure(1)
subplot(1,2,1)
plot(sort([StaticParams.aGrid,StaticParams.aGrid(1:end-1)]),Poly.pdf{2}(1,sort([1:length(Poly.pdf{2}),2:length(Poly.pdf{2})])),'r');
axis([StaticParams.a_min 10 0 0.62])
title('\xi^a_1');
xlabel('individual net worth a');
ylabel('probability');

subplot(1,2,2)
plot(sort([StaticParams.aGrid,StaticParams.aGrid(1:end-1)]),Poly.pdf{3}(1,sort([1:length(Poly.pdf{2}),2:length(Poly.pdf{2})])),'b');
axis([StaticParams.a_min 10 0 0.02])
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
axis([StaticParams.a_min 6 -7 7])
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
axis([StaticParams.a_min 6 -7 7])
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
        'Labels',{'PFI0','PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
% axis([0.5 8.5 -0.06 0.18])
subplot(1,2,2)
boxplot(a2,b2,...
        'Labels',{'PFI0','PFI1','PFI2','PFI3'},...,'PFI3'
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
        'Labels',{'PFI0','PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
axis([0.5 4.5 -4 -1])
subplot(1,2,2)
boxplot(a2,b2,...
        'Labels',{'PFI0','PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel B')   
axis([0.5 4.5 -4 -1])

a1=max(a1,[Res.stErr_PFI1_H;Res.stErr_PFI2_H;Res.stErr_PFI3_H;Res.stErr_PFI4_H]);
b1=[zeros(length(Res.stErr_PFI1_H),1);ones(length(Res.stErr_PFI2_H),1);...
    2*ones(length(Res.stErr_PFI3_H),1);3*ones(length(Res.stErr_PFI4_H),1)];
a2=max(a2,[Res.dynErr_PFI1_H;Res.dynErr_PFI2_H;Res.dynErr_PFI3_H;Res.dynErr_PFI4_H]);
b2=[zeros(length(Res.dynErr_PFI1_H),1);ones(length(Res.dynErr_PFI2_H),1);...
    2*ones(length(Res.dynErr_PFI3_H),1);3*ones(length(Res.dynErr_PFI4_H),1)];  

figure(5)
subplot(1,2,1)
boxplot(a1,b1,...
        'Labels',{'PFI0','PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel A')
axis([0.5 4.5 -2.5 0])
subplot(1,2,2)
boxplot(a2,b2,...
        'Labels',{'PFI0','PFI1','PFI2','PFI3'},...
        'Whisker',4,'OutlierSize',1,'Symbol','r.');
title('Panel B')   
axis([0.5 4.5 -2.5 0])



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
plot(StaticParams.aGrid,p*(Sol1.pdf-Sol0.pdf),'b')
plot(StaticParams.aGrid,p*(Sol2.pdf-Sol1.pdf),'r')
plot(StaticParams.aGrid,p*(Sol3.pdf-Sol2.pdf),'g')
plot([StaticParams.a_min StaticParams.a_max],StaticParams.criter_a*ones(1,2),'k',[StaticParams.a_min StaticParams.a_max],-StaticParams.criter_a*ones(1,2),'k')
hold off 
axis([StaticParams.a_min 30 -6e-4 6e-4])
legend('1^{st} order - 0^{th} order','2^{nd} order - 1^{st} order','3^{rd} order - 2^{nd} order','termination criterion','Location','southeast')
xlabel('individual net worth a');
ylabel('change in probability');
title('Panel A')
subplot(1,2,2)
hold on
plot(StaticParams.aGrid,p*(Sol1.pdf-Sol0.pdf),'b')
plot(StaticParams.aGrid,p*(Sol2.pdf-Sol1.pdf),'r')
plot(StaticParams.aGrid,p*(Sol3.pdf-Sol2.pdf),'g')
plot([StaticParams.a_min StaticParams.a_max],StaticParams.criter_a*ones(1,2),'k',[StaticParams.a_min StaticParams.a_max],-StaticParams.criter_a*ones(1,2),'k')
hold off 
axis([StaticParams.a_min 1.5 -1.8e-3 7e-4])
xlabel('individual net worth a');
ylabel('change in probability');
title('Panel B')


figure(7)
subplot(1,2,1)
hold on
plot(StaticParams.aGrid,p*Sol0.pdf,'r','LineWidth',1.5);
plot(StaticParams.aGrid,p*Sol1.pdf,'c','LineWidth',1.5);
plot(StaticParams.aGrid,p*Sol2.pdf,'b','LineWidth',1.5);
plot(StaticParams.aGrid,p*Sol3.pdf,'m--','LineWidth',1.5);
hold off
axis([StaticParams.a_min 15 0 0.009])
legend('PFI 0^{th} order','PFI 1^{st} order','PFI 2^{nd} order','PFI 3^{rd} order')
xlabel('individual net worth a');
ylabel('probability');
title('Panel A')
subplot(1,2,2)
hold on
x = sort([StaticParams.aGrid(1:end),StaticParams.aGrid(1:end-1)+1e-4]);
y = p*Sol0.pdf;
plot(x,y(sort([2:end,1:end])),'r','LineWidth',1.5);
y = p*Sol1.pdf;
plot(x,y(sort([2:end,1:end])),'c','LineWidth',1.5);
y = p*Sol2.pdf;
plot(x,y(sort([2:end,1:end])),'b','LineWidth',1.5);
y = p*Sol3.pdf;
plot(x,y(sort([2:end,1:end])),'m--','LineWidth',1.5);
hold off
axis([StaticParams.a_min 0.08 0.2 0.26])
xlabel('individual net worth a');
ylabel('probability');
title('Panel B')


%--------------------------------------------------------------------------
% Price Distribution
%--------------------------------------------------------------------------
zChange=[0;abs(Res.z_ag(2:end)-Res.z_ag(1:end-1))];
[a0,b0]=hist(squeeze(Res.q(end,logical((Res.z_ag==0).*zChange))),30);
[a1,b1]=hist(squeeze(Res.q(end,logical((Res.z_ag==0).*not(zChange)))),30);
[a2,b2]=hist(squeeze(Res.q(end,logical((Res.z_ag==1).*zChange))),30);
[a3,b3]=hist(squeeze(Res.q(end,logical((Res.z_ag==1).*not(zChange)))),30);

[c0,d0]=hist(squeeze(Res.r(end,logical((Res.z_ag==0).*zChange))),30);
[c1,d1]=hist(squeeze(Res.r(end,logical((Res.z_ag==0).*not(zChange)))),30);
[c2,d2]=hist(squeeze(Res.r(end,logical((Res.z_ag==1).*zChange))),30);
[c3,d3]=hist(squeeze(Res.r(end,logical((Res.z_ag==1).*not(zChange)))),30);

u=[0.025,((1+Res.r(end,2:end)).*Res.q(end,1:end-1)./Res.q(end,2:end)-1)];
[e0,f0]=hist((u(logical((Res.z_ag==0).*zChange))),30);
[e1,f1]=hist((u(logical((Res.z_ag==0).*not(zChange)))),30);
[e2,f2]=hist((u(logical((Res.z_ag==1).*zChange))),30);
[e3,f3]=hist((u(logical((Res.z_ag==1).*not(zChange)))),30);

figure(8)
subplot(1,3,1)
hold on
h0=plot([b0(1),sort([b0,b0(2:end)-1e-8]),b0(end)],[0,a0(sort([1:end,1:end-1])),0],'r');
h1=plot([b1(1),sort([b1,b1(2:end)-1e-8]),b1(end)],[0,a1(sort([1:end,1:end-1])),0],'m');
h2=plot([b2(1),sort([b2,b2(2:end)-1e-8]),b2(end)],[0,a2(sort([1:end,1:end-1])),0],'g');
h3=plot([b3(1),sort([b3,b3(2:end)-1e-8]),b3(end)],[0,a3(sort([1:end,1:end-1])),0],'c');
hold off
xlabel('house price');
ylabel('probability');
legend([h0,h1,h2,h3],{'first recession','repeated recession','first boom','repeated boom'},'Interpreter','latex','Location','North')
subplot(1,3,2)
hold on
plot([d0(1),sort([d0,d0(2:end)-1e-8]),d0(end)],[0,c0(sort([1:end,1:end-1])),0],'r');
plot([d1(1),sort([d1,d1(2:end)-1e-8]),d1(end)],[0,c1(sort([1:end,1:end-1])),0],'m');
plot([d2(1),sort([d2,d2(2:end)-1e-8]),d2(end)],[0,c2(sort([1:end,1:end-1])),0],'g');
plot([d3(1),sort([d3,d3(2:end)-1e-8]),d3(end)],[0,c3(sort([1:end,1:end-1])),0],'c');
hold off
xlabel('interest rate');
ylabel('probability');
subplot(1,3,3)
hold on
plot([f0(1),sort([f0,f0(2:end)-1e-8]),f0(end)],[0,e0(sort([1:end,1:end-1])),0],'r');
plot([f1(1),sort([f1,f1(2:end)-1e-8]),f1(end)],[0,e1(sort([1:end,1:end-1])),0],'m');
plot([f2(1),sort([f2,f2(2:end)-1e-8]),f2(end)],[0,e2(sort([1:end,1:end-1])),0],'g');
plot([f3(1),sort([f3,f3(2:end)-1e-8]),f3(end)],[0,e3(sort([1:end,1:end-1])),0],'c');
hold off
xlabel('user cost of housing');
ylabel('probability');


%--------------------------------------------------------------------------
% Policies
%--------------------------------------------------------------------------
b_prime = [getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.b(:,1,:)),Sol3.idx,Sol3.ergodic_weights(1,:));...
                 getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.b(:,2,:)),Sol3.idx,Sol3.ergodic_weights(1,:));...
                 getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.b(:,3,:)),Sol3.idx,Sol3.ergodic_weights(2,:));...
                 getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.b(:,4,:)),Sol3.idx,Sol3.ergodic_weights(2,:))];
h = max(0,[getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.h(:,1,:)),Sol3.idx,Sol3.ergodic_weights(1,:));...
                 getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.h(:,2,:)),Sol3.idx,Sol3.ergodic_weights(1,:));...
                 getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.h(:,3,:)),Sol3.idx,Sol3.ergodic_weights(2,:));...
                 getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.h(:,4,:)),Sol3.idx,Sol3.ergodic_weights(2,:))]);
c = max(0,[getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.c(:,1,:)),Sol3.idx,Sol3.ergodic_weights(1,:));...
                 getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.c(:,2,:)),Sol3.idx,Sol3.ergodic_weights(1,:));...
                 getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.c(:,3,:)),Sol3.idx,Sol3.ergodic_weights(2,:));...
                 getCApprox(Sol3.grid_N,Sol3.distrGrid,squeeze(Sol3.c(:,4,:)),Sol3.idx,Sol3.ergodic_weights(2,:))]);

figure(9)
subplot(2,2,1)
hold on
plot(StaticParams.aGrid_pol,c(1,:),'r');
plot(StaticParams.aGrid_pol,c(2,:),'m');
plot(StaticParams.aGrid_pol,c(3,:),'g');
plot(StaticParams.aGrid_pol,c(4,:),'c');
hold off
axis([StaticParams.a_min 15 0 1.8])
xlabel('individual net worth a');
ylabel('consumption policy');
legend('recession, low earner','recession, high earner','boom, low earner','boom, high earner','Location','southeast')
subplot(2,2,2)
hold on
plot(StaticParams.aGrid_pol,h(1,:),'r');
plot(StaticParams.aGrid_pol,h(2,:),'m');
plot(StaticParams.aGrid_pol,h(3,:),'g');
plot(StaticParams.aGrid_pol,h(4,:),'c');
hold off
axis([StaticParams.a_min 15 0 1.8])
xlabel('individual net worth a');
ylabel('housing policy');
subplot(2,2,3)
hold on
plot(StaticParams.aGrid_pol,b_prime(1,:),'r');
plot(StaticParams.aGrid_pol,b_prime(2,:),'m');
plot(StaticParams.aGrid_pol,b_prime(3,:),'g');
plot(StaticParams.aGrid_pol,b_prime(4,:),'c');
hold off
axis([StaticParams.a_min 15 -3.5 10.5])
xlabel('individual net worth a');
ylabel('bond policy');
subplot(2,2,4)
hold on
plot(StaticParams.aGrid_pol,b_prime(1,:),'r');
plot(StaticParams.aGrid_pol,b_prime(2,:),'m');
plot(StaticParams.aGrid_pol,b_prime(3,:),'g');
plot(StaticParams.aGrid_pol,b_prime(4,:),'c');
hold off
axis([StaticParams.a_min 1 -3.5 0.5])
xlabel('individual net worth a');
ylabel('bond policy');

%--------------------------------------------------------------------------
% Compute Times
%--------------------------------------------------------------------------
Sol0 = load(strcat(folder,'/Sol0_PFI.mat'));
Sol1 = load(strcat(folder,'/Sol1_PFI.mat'));
Sol2 = load(strcat(folder,'/Sol2_PFI.mat'));
Sol3 = load(strcat(folder,'/Sol3_PFI.mat'));
preSol = load(strcat(folder,'/preSol_PFI.mat'));

T_PFI = preSol.time_b+preSol.time_g+preSol.time_distr;
T_PFI = T_PFI+cumsum([Sol0.time_solve;Sol1.time_solve;Sol2.time_solve;Sol3.time_solve]);
T_PFI = T_PFI./60;
T_PFI_wc = preSol.time_b_wallclock+preSol.time_g_wallclock+preSol.time_distr_wallclock;
T_PFI_wc = T_PFI_wc+cumsum([Sol0.time_solve_wallclock;Sol1.time_solve_wallclock;Sol2.time_solve_wallclock;Sol3.time_solve_wallclock]);
T_PFI_wc = T_PFI_wc./60;

%--------------------------------------------------------------------------
% Ergodic Gini
%--------------------------------------------------------------------------

pdf = 0.5.*[sum(Sol0.pdf,1);sum(Sol1.pdf,1);sum(Sol2.pdf,1);sum(Sol3.pdf,1)];
ergodic_gini_PFI = [sum(sum((pdf(1,:)'*pdf(1,:)).*abs(repmat(StaticParams.aGrid',[1,size(pdf,2)])-repmat(StaticParams.aGrid,[size(pdf,2),1])),1),2)...
            /(2*StaticParams.aGrid*pdf(1,:)');...
            sum(sum((pdf(2,:)'*pdf(2,:)).*abs(repmat(StaticParams.aGrid',[1,size(pdf,2)])-repmat(StaticParams.aGrid,[size(pdf,2),1])),1),2)...
            /(2*StaticParams.aGrid*pdf(2,:)');...
            sum(sum((pdf(3,:)'*pdf(3,:)).*abs(repmat(StaticParams.aGrid',[1,size(pdf,2)])-repmat(StaticParams.aGrid,[size(pdf,2),1])),1),2)...
            /(2*StaticParams.aGrid*pdf(3,:)');...
            sum(sum((pdf(4,:)'*pdf(4,:)).*abs(repmat(StaticParams.aGrid',[1,size(pdf,2)])-repmat(StaticParams.aGrid,[size(pdf,2),1])),1),2)...
            /(2*StaticParams.aGrid*pdf(4,:)')];
