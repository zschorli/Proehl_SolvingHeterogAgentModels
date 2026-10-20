% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This function solves the model with aggregate risk.
%__________________________________________________________________________
function solveModel_AggShock(StaticParams,folder,order)

mu = 1e-2;
SymFuns = buildSymFuns(StaticParams, mu);

% Vary the different truncation orders of the polynomial chaos expansion.
% Start with order 1 as order zero is fixed due to the zero net supply
% condition.
for M=0:order
    
%--------------------------------------------------------------------------
% Set the grid for polynomial chaos coefficients
%--------------------------------------------------------------------------
Poly =load(strcat(folder,'/Poly',num2str(M),'_PFI.mat'));

Sol = struct;
grid_N1 = 7;
grid_N2 = 3^(Poly.phi_N>0);
grid_N3 = 3^(Poly.phi_N>1);
grid_N4 = 2^(Poly.phi_N>2); 

Sol.distrGrid_min = [2.9, 1.1, 7e-2, -3e-2];
Sol.distrGrid_max = [3.3, 1.3, 2e-1, -1e-2];

Sol.grid_N = [grid_N1, grid_N2, grid_N3, grid_N4];
G = zeros(4,grid_N4,grid_N3,grid_N2,grid_N1);
[G(4,:,:,:,:,:),G(3,:,:,:,:,:),G(2,:,:,:,:,:),G(1,:,:,:,:,:)] = ...
ndgrid(1:Sol.grid_N(4),1:Sol.grid_N(3),1:Sol.grid_N(2),1:Sol.grid_N(1));
grid_N = grid_N1*grid_N2*grid_N3*grid_N4;

Sol.distrGrid = zeros(max([2,grid_N1,grid_N2,grid_N3,grid_N4]),Poly.phi_N+1);
Sol.idx = zeros(grid_N,Poly.phi_N+1);
for i=1:Poly.phi_N+1
    Sol.distrGrid(1:Sol.grid_N(i),i) = linspace(Sol.distrGrid_min(i),Sol.distrGrid_max(i),Sol.grid_N(i))';
    Sol.idx(:,i) = reshape(squeeze(G(i,:,:,:,:,:)),[],1);
end

%--------------------------------------------------------------------------
% Initialization 
%--------------------------------------------------------------------------
Sol.r = zeros(grid_N,2); 
% Note that q equals the mean of the savings distribution which equals the
% zeroth-order projection coefficient Sol.q = Sol.distrGrid(Sol.idx(:,1),1)
Sol.wealth = zeros(grid_N,4,length(StaticParams.aGrid_pol));
Sol.c = zeros(grid_N,4,length(StaticParams.aGrid_pol));
Sol.h = zeros(grid_N,4,length(StaticParams.aGrid_pol));
Sol.a_prime = zeros(grid_N,4,length(StaticParams.aGrid_pol));
Sol.diff = 1e8*ones(1,3);
grid = zeros(grid_N,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:)),length(Poly.xi{3}(1,:)));

if M<=0
    preSol =load(strcat(folder,'/preSol_PFI.mat'));
    [Yz,Ya] = ndgrid(1:2,StaticParams.aGrid_pol);
    c_g_approx = @(r,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,preSol.c_g,Yz,Ya,(1+r).*ones(size(Yz)),q.*ones(size(Yz)));
    c_b_approx = @(r,q) interpn(1:2,StaticParams.aGrid_pol,StaticParams.RGrid,StaticParams.qGrid,preSol.c_b,Yz,Ya,(1+r).*ones(size(Yz)),q.*ones(size(Yz)));
    c_ss = [c_b_approx(preSol.agR_stationary_b,preSol.agA_stationary_b);...
            c_g_approx(preSol.agR_stationary_g,preSol.agA_stationary_g)];
else
    preSol = load(strcat(folder,'/Sol',num2str(M-1),'_PFI.mat')); 
end

for m=1:grid_N
    d = num2cell(Sol.idx(m,:)); 
    weights = diag(Sol.distrGrid([d{:}],:))';
    grid(m,:,:,:,:) = squeeze(sum(repmat(weights',[1,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:)),length(Poly.xi{3}(1,:))]).*Poly.total,1));
    if M<=0     
        stationaryUq = [preSol.agR_stationary_b,preSol.agR_stationary_g]...
            .*[preSol.agA_stationary_b,preSol.agA_stationary_g];
        Sol.q_pr(m,:) = repmat(Sol.distrGrid(Sol.idx(m,1),1),[1,2]);
        Sol.r(m,:) = ((stationaryUq+Sol.q_pr(m,:))./repmat(Sol.distrGrid(Sol.idx(m,1),1),[1,2])-1);
        Sol.wealth(m,:,:) = StaticParams.wage...
                        +repmat(1+Sol.r(m,[1,1,2,2])',[1,length(StaticParams.aGrid_pol)])... 
                        .*repmat(StaticParams.aGrid_pol,[4,1]);
        Sol.c(m,:,:) = c_ss;
    end
end
if M<=0
    r_large = repmat(Sol.r(:,[1,1,2,2]),[1,1,size(Sol.c,3)]);
    q_pr_large = repmat(Sol.q_pr(:,[1,1,2,2]),[1,1,size(Sol.c,3)]);
    q = repmat(Sol.distrGrid(Sol.idx(:,1),1),[1,4,size(Sol.c,3)]);
    a_large = repmat(permute(StaticParams.aGrid_pol,[3,1,2]),[size(Sol.c,1),4,1]);
    uq = max(eps,(1+r_large).*q-q_pr_large);
    h_bound = (1/(1-StaticParams.theta).*a_large./q);
    Sol.h = min((StaticParams.eta.*(uq).^(-StaticParams.eps).*Sol.c),h_bound);
    Sol.a_prime = max(0,Sol.wealth-Sol.c-(uq).*Sol.h);
    Sol.a_prime = Sol.a_prime./max(max(max(Sol.a_prime))).*(0.99.*StaticParams.a_max+0.*max(max(max(Sol.a_prime))));
    Sol.c = Sol.wealth-Sol.a_prime-(uq).*Sol.h;
else
    aux = reshape(repmat(1:size(preSol.a_prime,1),[Sol.grid_N(M+1),1]),[],1);
    Sol.q_pr = preSol.q_pr(aux,:);
    Sol.r = preSol.r(aux,:);
    Sol.wealth = preSol.wealth(aux,:,:);
    Sol.a_prime = preSol.a_prime(aux,:,:);
    Sol.c = preSol.c(aux,:,:);
    Sol.h = preSol.h(aux,:,:);
end
clear a_prime_g_approx a_prime_b_approx c_g_approx c_b_approx...
      y_g_approx y_b_approx d i j m unemplGrid emplGrid;

%--------------------------------------------------------------------------
% Iteration 
%--------------------------------------------------------------------------
tic;
tstart = cputime;
[Sol.c,Sol.h,Sol.a_prime,Sol.b,Sol.wealth,Sol.r,Sol.q_pr,Sol.diff]...
= OptPol_AggShock_PFI(Sol.c,Sol.h,Sol.a_prime,Sol.wealth,Sol.r,Sol.q_pr,...
  grid,StaticParams,Poly,Sol.grid_N,Sol.distrGrid,Sol.idx,SymFuns);
Sol.time_solve_wallclock = toc;
Sol.time_solve = cputime-tstart;

save(strcat(folder,'/Sol',num2str(M),'_PFI.mat'),'-struct','Sol'); 
end
end