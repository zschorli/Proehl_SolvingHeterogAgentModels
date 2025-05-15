% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function solves the model with aggregate risk.
%__________________________________________________________________________
function solveModel_AggShock(StaticParams,case_nr,folder,order)

% Vary the different truncation orders of the polynomial chaos expansion.
for M=0:order
    
%--------------------------------------------------------------------------
% Set the grid for polynomial chaos coefficients
%--------------------------------------------------------------------------
Poly =load(strcat(folder,'/Poly',num2str(M),'_PFI.mat'));

Sol = struct;
grid_N1 = 6;
grid_N2 = 11^(Poly.phi_N>0);
grid_N3 = 9^(Poly.phi_N>1);
grid_N4 = 11^(Poly.phi_N>2);
grid_N5 = 2^(Poly.phi_N>3);
if StaticParams.rho==0
    Sol.distrGrid_min = [10, 1,  -5e-3,-2e-3, -1e-6];
    Sol.distrGrid_max = [14, 1.6, 5e-2, 3e-3,  1e-6];
elseif StaticParams.k_min>0
    Sol.distrGrid_min = [20,  1e-2, -2e-2, -1e-4, -1e-6];
    Sol.distrGrid_max = [23.5, 0.4,  1e-4,  8e-4,  1e-6];
else
    Sol.distrGrid_min = [10.5, 1, 1e-3, -2e-3, -1e-6];
    Sol.distrGrid_max = [13, 1.2, 3e-2,  3e-4,  1e-6];
end
Sol.grid_N = [grid_N1, grid_N2, grid_N3, grid_N4, grid_N5];
G = zeros(5,grid_N5,grid_N4,grid_N3,grid_N2,grid_N1);
[G(5,:,:,:,:,:),G(4,:,:,:,:,:),G(3,:,:,:,:,:),G(2,:,:,:,:,:),G(1,:,:,:,:,:)] = ...
ndgrid(1:Sol.grid_N(5),1:Sol.grid_N(4),1:Sol.grid_N(3),1:Sol.grid_N(2),1:Sol.grid_N(1));
grid_N = grid_N1*grid_N2*grid_N3*grid_N4*grid_N5;

Sol.distrGrid = zeros(max(Sol.grid_N),Poly.phi_N+1);
Sol.idx = zeros(grid_N,Poly.phi_N+1);
for i=1:Poly.phi_N+1
    Sol.distrGrid(1:Sol.grid_N(i),i) = linspace(Sol.distrGrid_min(i),Sol.distrGrid_max(i),Sol.grid_N(i))';
    Sol.idx(:,i) = reshape(squeeze(G(i,:,:,:,:,:)),[],1);
end
 
clear i;

%--------------------------------------------------------------------------
% Initialization 
%--------------------------------------------------------------------------
Sol.agK = zeros(grid_N,1);
Sol.agK1 = zeros(grid_N,2);
Sol.agK2 = zeros(grid_N,2);
Sol.agL = zeros(grid_N,2);
Sol.agLtax = zeros(grid_N,2);
Sol.rate = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.wage = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.wealth = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.c = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.l = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.k_prime = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.diff = 1e8*ones(1,2);

grid = zeros(grid_N,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:)));

if M==0
    switch case_nr
    case {1,5,6,13,14,19,20,21}
        count = 1;
    case {2,3,4,7}
        count = case_nr;
    otherwise
        count = 7;
    end
    preSol =load(strcat('res_case',num2str(count),'/preSol_PFI.mat'));
    if numel(StaticParams.LtaxGrid)>1
        c_g_approx = @(z_idx,k,K,L,Ltax) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,squeeze(preSol.c_g(z_idx,:,:,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)),min(StaticParams.LGrid(end),max(StaticParams.LGrid(1),L)),min(StaticParams.LtaxGrid(end),max(StaticParams.LtaxGrid(1),Ltax)));
        c_b_approx = @(z_idx,k,K,L,Ltax) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,squeeze(preSol.c_b(z_idx,:,:,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)),min(StaticParams.LGrid(end),max(StaticParams.LGrid(1),L)),min(StaticParams.LtaxGrid(end),max(StaticParams.LtaxGrid(1),Ltax)));
        l_g_approx = @(z_idx,k,K,L,Ltax) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,squeeze(preSol.l_g(z_idx,:,:,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)),min(StaticParams.LGrid(end),max(StaticParams.LGrid(1),L)),min(StaticParams.LtaxGrid(end),max(StaticParams.LtaxGrid(1),Ltax)));
        l_b_approx = @(z_idx,k,K,L,Ltax) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,StaticParams.LtaxGrid,squeeze(preSol.l_b(z_idx,:,:,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)),min(StaticParams.LGrid(end),max(StaticParams.LGrid(1),L)),min(StaticParams.LtaxGrid(end),max(StaticParams.LtaxGrid(1),Ltax)));
    else
        c_g_approx = @(z_idx,k,K,L,Ltax) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,squeeze(preSol.c_g(z_idx,:,:,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)),min(StaticParams.LGrid(end),max(StaticParams.LGrid(1),L)));
        c_b_approx = @(z_idx,k,K,L,Ltax) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,squeeze(preSol.c_b(z_idx,:,:,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)),min(StaticParams.LGrid(end),max(StaticParams.LGrid(1),L)));
        l_g_approx = @(z_idx,k,K,L,Ltax) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,squeeze(preSol.l_g(z_idx,:,:,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)),min(StaticParams.LGrid(end),max(StaticParams.LGrid(1),L)));
        l_b_approx = @(z_idx,k,K,L,Ltax) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,StaticParams.LGrid,squeeze(preSol.l_b(z_idx,:,:,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)),min(StaticParams.LGrid(end),max(StaticParams.LGrid(1),L)));
    end
    c_ss = [c_b_approx(1,StaticParams.kGrid_pol,preSol.agK_stationary_b,preSol.agL_stationary_b,preSol.agLtaxRatio_stationary_b)';...
            c_b_approx(2,StaticParams.kGrid_pol,preSol.agK_stationary_b,preSol.agL_stationary_b,preSol.agLtaxRatio_stationary_b)';...
            c_g_approx(1,StaticParams.kGrid_pol,preSol.agK_stationary_g,preSol.agL_stationary_g,preSol.agLtaxRatio_stationary_g)';...
            c_g_approx(2,StaticParams.kGrid_pol,preSol.agK_stationary_g,preSol.agL_stationary_g,preSol.agLtaxRatio_stationary_g)'];
    l_ss = [l_b_approx(1,StaticParams.kGrid_pol,preSol.agK_stationary_b,preSol.agL_stationary_b,preSol.agLtaxRatio_stationary_b)';...
            l_b_approx(2,StaticParams.kGrid_pol,preSol.agK_stationary_b,preSol.agL_stationary_b,preSol.agLtaxRatio_stationary_b)';...
            l_g_approx(1,StaticParams.kGrid_pol,preSol.agK_stationary_g,preSol.agL_stationary_g,preSol.agLtaxRatio_stationary_g)';...
            l_g_approx(2,StaticParams.kGrid_pol,preSol.agK_stationary_g,preSol.agL_stationary_g,preSol.agLtaxRatio_stationary_g)'];
else
    preSol = load(strcat(folder,'/Sol',num2str(M-1),'_PFI.mat')); 
    aux = reshape(repmat(1:size(preSol.k_prime,1),[Sol.grid_N(M+1),1]),[],1);
    Sol.k_prime = preSol.k_prime(aux,:,:);
    Sol.c = preSol.c(aux,:,:);
    Sol.l = preSol.l(aux,:,:);
end

for m=1:grid_N
    d = num2cell(Sol.idx(m,:)); 
    weights = diag(Sol.distrGrid([d{:}],:))';
    grid(m,:,:) = squeeze(sum(repmat(weights',[1,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:))]).*Poly.total,1));
           
    Sol.agK(m) = sum(sum(Poly.total_pdf.*squeeze(grid(m,:,:)),2),1);    
	Sol.agK1(m,:) = [sum(sum(Poly.total_pdf(1:2,:).*squeeze(grid(m,1:2,:)),2),1)/StaticParams.ur(1);...
                     sum(squeeze(Poly.total_pdf(1,:)').*squeeze(grid(m,1,:)),1)/StaticParams.ur(2)];
	Sol.agK2(m,:) = [sum(squeeze(Poly.total_pdf(3,:)').*squeeze(grid(m,3,:)),1)/StaticParams.er(1);...
                     sum(sum(Poly.total_pdf(2:3,:).*squeeze(grid(m,2:3,:)),2),1)/StaticParams.er(2)];
    if M==0
        Sol.l(m,:,:) = l_ss;
    end
    l_prime_large = interp1(StaticParams.kGrid_pol,squeeze(Sol.l(m,:,:))',StaticParams.kGrid)';   
    [~,~,~,Sol.agL(m,:)] = calcAggregates(l_prime_large,squeeze(grid(m,:,:)),StaticParams,Poly);
    [~,~,~,Sol.agLtax(m,:)] = calcAggregates(l_prime_large.^(1-StaticParams.rho),squeeze(grid(m,:,:)),StaticParams,Poly);
    
    % interest rate 
    Sol.rate(m,:,:) = [repmat(StaticParams.rate(0,Sol.agK(m),Sol.agL(m,1)),[1,numel(StaticParams.kGrid_pol)]);...
                       repmat(StaticParams.rate(1,Sol.agK(m),Sol.agL(m,2)),[1,numel(StaticParams.kGrid_pol)])];
    % wage 
    Sol.wage(m,:,:) = [repmat(StaticParams.wage(0,Sol.agK(m),Sol.agL(m,1)),[1,numel(StaticParams.kGrid_pol)]);...
                       repmat(StaticParams.wage(1,Sol.agK(m),Sol.agL(m,2)),[1,numel(StaticParams.kGrid_pol)])];   
    % income
    Sol.wealth(m,:,:) = [StaticParams.wealth(0,Sol.agK(m),Sol.agL(m,1),Sol.agL(m,1)./Sol.agLtax(m,1),StaticParams.kGrid_pol,squeeze(Sol.l(m,1:2,:)));...
                         StaticParams.wealth(1,Sol.agK(m),Sol.agL(m,2),Sol.agL(m,2)./Sol.agLtax(m,2),StaticParams.kGrid_pol,squeeze(Sol.l(m,3:4,:)))]; 
    if M==0
        Sol.c(m,:,:) = min(squeeze(Sol.wealth(m,:,:)),c_ss);
        Sol.k_prime(m,:,:) = squeeze(Sol.wealth(m,:,:)-Sol.c(m,:,:));
    end
end

clear k_prime_g_approx k_prime_b_approx c_g_approx c_b_approx...
      y_g_approx y_b_approx d i j m unemplGrid emplGrid;

%--------------------------------------------------------------------------
% Iteration 
%--------------------------------------------------------------------------
tic;
[Sol.c,Sol.l,Sol.k_prime,Sol.agL,Sol.agLtax,Sol.diff]...
= OptPol_AggShock_PFI(Sol.c,Sol.l,Sol.k_prime,Sol.agL,Sol.agLtax,Sol.agK,grid,StaticParams,Poly,Sol.grid_N,Sol.distrGrid,Sol.idx);
Sol.time_solve = toc;

for m=1:grid_N
    % interest rate 
    Sol.rate(m,:,:) = [repmat(StaticParams.rate(0,Sol.agK(m),Sol.agL(m,1)),[1,numel(StaticParams.kGrid_pol)]);...
                       repmat(StaticParams.rate(1,Sol.agK(m),Sol.agL(m,2)),[1,numel(StaticParams.kGrid_pol)])];
    % wage 
    Sol.wage(m,:,:) = [repmat(StaticParams.wage(0,Sol.agK(m),Sol.agL(m,1)),[1,numel(StaticParams.kGrid_pol)]);...
                       repmat(StaticParams.wage(1,Sol.agK(m),Sol.agL(m,2)),[1,numel(StaticParams.kGrid_pol)])];   
    % income
    Sol.wealth(m,:,:) = [StaticParams.wealth(0,Sol.agK(m),Sol.agL(m,1),Sol.agL(m,1)./Sol.agLtax(m,1),StaticParams.kGrid_pol,squeeze(Sol.l(m,1:2,:)));...
                         StaticParams.wealth(1,Sol.agK(m),Sol.agL(m,2),Sol.agL(m,2)./Sol.agLtax(m,2),StaticParams.kGrid_pol,squeeze(Sol.l(m,3:4,:)))]; 
end
Sol.k_prime = Sol.k_prime.*(Sol.k_prime>eps);
save(strcat(folder,'/Sol',num2str(M),'_PFI.mat'),'-struct','Sol'); 
end
end