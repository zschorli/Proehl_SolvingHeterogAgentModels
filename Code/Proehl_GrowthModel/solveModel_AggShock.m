% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This function solves the model with aggregate risk.
%__________________________________________________________________________
function solveModel_AggShock(StaticParams,case_nr,folder,order)

% Use either the Proximal point algorithm (method = 0) or policy function
% iteration (method = 1)
for method = 1:-1:0
StaticParams.methodEqPFI = method;  

% Vary the different truncation orders of the polynomial chaos expansion.
% Start with order 1 as order zero is fixed due to the zero net supply
% condition.
for M=0:order
    
%--------------------------------------------------------------------------
% Set the grid for polynomial chaos coefficients
%--------------------------------------------------------------------------
if StaticParams.methodEqPFI == 0
	Poly =load(strcat(folder,'/Poly',num2str(M),'_PPA.mat'));
else
	Poly =load(strcat(folder,'/Poly',num2str(M),'_PFI.mat'));
end
grid_N1 = 4;
grid_N2 = 3^(Poly.phi_N>0);
grid_N3 = 2^(Poly.phi_N>1);
grid_N4 = 2^(Poly.phi_N>2);
grid_N = grid_N1*grid_N2*grid_N3*grid_N4^(Poly.phi_N-2);

Sol = struct;
switch case_nr
    case {2,4}
        Sol.distrGrid_min = [36.5, 0.99, -1e-4, -1e-6, -5e-08];
        Sol.distrGrid_max = [43  , 1.2,   5e-3,  5e-6,  5e-08];
    case {5,13}
        Sol.distrGrid_min = [3, 0.86, -1.5e-2, -4e-3];
        Sol.distrGrid_max = [8, 1,     2.5e-2,  5e-5];
    case {6,14}
        Sol.distrGrid_min = [52, 0.96,-1e-3, -8e-5];
        Sol.distrGrid_max = [68, 1.3,  4e-4,  1e-5];
    case {9,17}
        Sol.distrGrid_min = [3.5, 0.7, -1e-2, -7e-3];
        Sol.distrGrid_max = [8,   1,    5e-2,  5e-5];
    case {10,18}
        Sol.distrGrid_min = [95,  0.96, -1e-3, -8e-5];
        Sol.distrGrid_max = [115, 1.3,   4e-4,  1e-5];
    case {11,19}
        Sol.distrGrid_min = [80,  0.95, -1e-3, -8e-5];
        Sol.distrGrid_max = [100, 1.3,   4e-4,  1e-5];
    case {12,20}
        Sol.distrGrid_min = [5,  0.83, -1.5e-2, -6e-3];
        Sol.distrGrid_max = [12, 1,     3e-2,    5e-5];
    otherwise
        Sol.distrGrid_min = [36.5, 0.99, -5e-4, -1e-6];
        Sol.distrGrid_max = [43  , 1.1 ,  5e-4,  5e-5];
end
Sol.distrGrid = zeros(grid_N1,Poly.phi_N+1);
Sol.distrGrid(1:grid_N1,1) = linspace(Sol.distrGrid_min(1),Sol.distrGrid_max(1),grid_N1)';  
if Poly.phi_N>0
    Sol.distrGrid(1:grid_N2,2) = linspace(Sol.distrGrid_min(2),Sol.distrGrid_max(2),grid_N2)';    
end
if Poly.phi_N>1
    Sol.distrGrid(1:grid_N3,3) = linspace(Sol.distrGrid_min(3),Sol.distrGrid_max(3),grid_N3)';
end
if Poly.phi_N>2
    for i=4:Poly.phi_N+1
        Sol.distrGrid(1:grid_N4,i) = linspace(Sol.distrGrid_min(i),Sol.distrGrid_max(i),grid_N4)';
    end
end

Sol.idx = zeros(grid_N,Poly.phi_N+1);
for j=0:grid_N/grid_N1:grid_N/grid_N1*(grid_N1-1)
    Sol.idx(j+1:j+grid_N/grid_N1,1) = 1+mod(j/(grid_N/grid_N1),grid_N1); 
end
if Poly.phi_N>0
    for j=0:grid_N/grid_N1/grid_N2:grid_N/grid_N1/grid_N2*(grid_N1*grid_N2-1)
        Sol.idx(j+1:j+grid_N/grid_N1/grid_N2,2) = 1+mod(j/(grid_N/grid_N1/grid_N2),grid_N2); 
    end
end
if Poly.phi_N>1
    for j=0:grid_N/grid_N1/grid_N2/grid_N3:grid_N/grid_N1/grid_N2/grid_N3*(grid_N1*grid_N2*grid_N3-1)
        Sol.idx(j+1:j+grid_N/grid_N1/grid_N2/grid_N3,3) = 1+mod(j/(grid_N/grid_N1/grid_N2/grid_N3),grid_N3); 
    end
end
if Poly.phi_N>2
    for i=1:(Poly.phi_N-2)
        for j=0:grid_N/grid_N1/grid_N2/grid_N3/grid_N4^i:grid_N/grid_N1/grid_N2/grid_N3/grid_N4^i*(grid_N4^i*grid_N1*grid_N2*grid_N3-1)
            Sol.idx(j+1:j+grid_N/grid_N1/grid_N2/grid_N3/grid_N4^i,3+i) = 1+mod(j/(grid_N/grid_N1/grid_N2/grid_N3/grid_N4^i),grid_N4); 
        end
    end
end
 
clear cdf_proj idx_proj kGrid kGrid i j;

%--------------------------------------------------------------------------
% Initialization 
%--------------------------------------------------------------------------
Sol.agK = zeros(grid_N,1);
Sol.agK1 = zeros(grid_N,2);
Sol.agK2 = zeros(grid_N,2);
Sol.rate = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.wage = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.wealth = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.c = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.k_prime = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.diff = 1e8*ones(1,2);

c_new = zeros(grid_N,4,length(StaticParams.kGrid_pol));
k_prime_new = zeros(grid_N,4,length(StaticParams.kGrid_pol));

if StaticParams.methodEqPFI == 0
	Sol.y = zeros(grid_N,4,length(StaticParams.kGrid_pol));
	Sol.sol_y = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	Sol.sol_nu = zeros(grid_N,4*length(StaticParams.kGrid_pol));
    Sol.sol_y_1 = zeros(grid_N,4,length(StaticParams.kGrid_pol));
	Sol.sol_y_2 = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	Sol.sol_nu_2 = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	Sol.A = 10*ones(1,grid_N);
	Sol.flag = ones(1,grid_N);

	y_new = zeros(grid_N,4,length(StaticParams.kGrid_pol));
	sol_y_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	sol_nu_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	sol_y_2_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	sol_nu_2_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	A_new = Sol.A;
	flag = zeros(grid_N,1);
end

grid = zeros(grid_N,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:)));%,length(Poly.xi{3}(1,:)));
diffs = 1e8*ones(grid_N,2);

if StaticParams.methodEqPFI == 0
	preSol =load(strcat(folder,'/preSol_PPA.mat'));
else
	preSol =load(strcat(folder,'/preSol_PFI.mat'));
end
c_g_approx = @(z_idx,k,K) interpn(StaticParams.kGrid_pol,preSol.KGrid,squeeze(preSol.c_g(z_idx,:,:)),k,min(preSol.KGrid(end),max(preSol.KGrid(1),K)));
c_b_approx = @(z_idx,k,K) interpn(StaticParams.kGrid_pol,preSol.KGrid,squeeze(preSol.c_b(z_idx,:,:)),k,min(preSol.KGrid(end),max(preSol.KGrid(1),K)));
K_ss = [(2*StaticParams.p(1:2)'*preSol.pdf_cond_b*StaticParams.kGrid');...
        (2*StaticParams.p(3:4)'*preSol.pdf_cond_g*StaticParams.kGrid')];
c_ss = sort([c_b_approx(1,StaticParams.kGrid_pol,K_ss(1))';...
             c_b_approx(2,StaticParams.kGrid_pol,K_ss(1))';...
             c_g_approx(1,StaticParams.kGrid_pol,K_ss(2))';...
             c_g_approx(2,StaticParams.kGrid_pol,K_ss(2))'],2);
for m=1:grid_N
    d = num2cell(Sol.idx(m,:)); 
    weights = diag(Sol.distrGrid([d{:}],:))';
    grid(m,:,:) = squeeze(sum(repmat(weights',[1,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:))]).*Poly.total,1));%,length(Poly.xi{3}(1,:))
           
    Sol.agK(m) = sum(sum(Poly.total_pdf.*squeeze(grid(m,:,:)),2),1);    
	Sol.agK1(m,:) = [sum(sum(Poly.total_pdf(1:2,:).*squeeze(grid(m,1:2,:)),2),1)/StaticParams.ur(1);...
                     sum(squeeze(Poly.total_pdf(1,:)').*squeeze(grid(m,1,:)),1)/StaticParams.ur(2)];
	Sol.agK2(m,:) = [sum(squeeze(Poly.total_pdf(3,:)').*squeeze(grid(m,3,:)),1)/StaticParams.er(1);...
                     sum(sum(Poly.total_pdf(2:3,:).*squeeze(grid(m,2:3,:)),2),1)/StaticParams.er(2)];
    
    % interest rate 
    Sol.rate(m,:,:) = StaticParams.alpha.*...
                        repmat([(1-StaticParams.delta_a)*ones(2,1);(1+StaticParams.delta_a)*ones(2,1)],[1,length(StaticParams.kGrid_pol)])...
                       .*([repmat(Sol.agK(m,1),[2,length(StaticParams.kGrid_pol)]);repmat(Sol.agK(m,1),[2,length(StaticParams.kGrid_pol)])]...
                       ./(StaticParams.l_bar.*repmat([StaticParams.er(1)*ones(2,1);StaticParams.er(2)*ones(2,1)],[1,length(StaticParams.kGrid_pol)])))...
                       .^(StaticParams.alpha-1);   
    % wage 
    Sol.wage(m,:,:) = (1-StaticParams.alpha).*...
                        repmat([(1-StaticParams.delta_a)*ones(2,1);(1+StaticParams.delta_a)*ones(2,1)],[1,length(StaticParams.kGrid_pol)])...
                       .*([repmat(Sol.agK(m,1),[2,length(StaticParams.kGrid_pol)]);repmat(Sol.agK(m,1),[2,length(StaticParams.kGrid_pol)])]...
                       ./(StaticParams.l_bar.*repmat([StaticParams.er(1)*ones(2,1);StaticParams.er(2)*ones(2,1)],[1,length(StaticParams.kGrid_pol)])))...
                       .^(StaticParams.alpha);   
    % income
    Sol.wealth(m,:,:) = ((StaticParams.l_bar-...
                          repmat([StaticParams.mu.*StaticParams.ur(1)/StaticParams.er(1)*ones(2,1);...
                                  StaticParams.mu.*StaticParams.ur(2)/StaticParams.er(2)*ones(2,1)],[1,length(StaticParams.kGrid_pol)]))...
                          .*repmat([0;1;0;1],[1,length(StaticParams.kGrid_pol)])...
                          +StaticParams.mu.*(1-repmat([0;1;0;1],[1,length(StaticParams.kGrid_pol)])))...
                          .*squeeze(Sol.wage(m,:,:))...
                          +(1-StaticParams.delta+squeeze(Sol.rate(m,:,:)))...
                          .*repmat(StaticParams.kGrid_pol,[4,1]);
    
    Sol.c(m,:,:) = min(squeeze(Sol.wealth(m,:,:)),c_ss);
    Sol.k_prime(m,:,:) = squeeze(Sol.wealth(m,:,:)-Sol.c(m,:,:));
    if StaticParams.methodEqPFI == 0
		Sol.y(m,:,:) = 0;
		Sol.sol_y(m,:) = reshape(squeeze(Sol.k_prime(m,:,:))',[4*length(StaticParams.kGrid_pol),1]);
		Sol.sol_nu(m,:) = Sol.sol_y(m,:);
		Sol.sol_y_1(m,:,:) = Sol.c(m,:,:);
		Sol.sol_y_2(m,:) = reshape(squeeze(Sol.y(m,:,:))',[4*length(StaticParams.kGrid_pol),1]);
		Sol.sol_nu_2(m,:) = Sol.sol_y_2(m,:);
    end
end

clear k_prime_g_approx k_prime_b_approx c_g_approx c_b_approx...
      y_g_approx y_b_approx d i j m unemplGrid emplGrid;

iter = 1;
%--------------------------------------------------------------------------
% Iteration 
%--------------------------------------------------------------------------
if StaticParams.methodEqPFI == 1
	tic;
	[Sol.c,Sol.k_prime,Sol.diff]...
    = OptPol_AggShock_PFI(Sol.c,Sol.k_prime,grid,StaticParams,Poly,Sol.distrGrid,Sol.idx,Sol.wealth);
    Sol.time_solve = toc;
else
	agK_prime_new_var = zeros(grid_N,2,2,length(StaticParams.kGrid_pol));
	c_pr_var = zeros(grid_N,4,4,length(StaticParams.kGrid_pol));
    
	tic;
	while max(Sol.diff(end,:)) > StaticParams.criter_k && sum(Sol.flag(end,:)>0)>0
        parfor m=1:grid_N  %
            k_prime = squeeze(Sol.k_prime(m,:,:));
			y = squeeze(Sol.y(m,:,:));
			sol_y = Sol.sol_y(m,:)';
			sol_y_2 = Sol.sol_y_2(m,:)';
			sol_nu = Sol.sol_nu(m,:)';
			sol_nu_2 = Sol.sol_nu_2(m,:)';
			A = Sol.A(m);
            
            % compute the derivatives of aggregate variables w.r.t. changes
            % in individual capital (needed for the Jacobian in the
            % nonlinear solver which solves the Euler equation)
            condition = iter<=3 ...
                        || sum(sum(sum(c_pr_var(m,:,:,:),2),3),4)==0 ...
                        || floor(max(sum(Sol.diff(3:end,:),1))/norm(3*Sol.diff(1,:)))...
                           >floor(max(sum(Sol.diff(3:end-1,:),1))/norm(3*Sol.diff(1,:)));
            if condition
				[agK_prime_new_var(m,:,:,:),c_pr_var(m,:,:,:)] = calcJacobianProxies(reshape(sol_y,size(k_prime'))',squeeze(grid(m,:,:,:)),StaticParams,Sol,Poly);
            end 
            
            % update of the proximal point algorithm
			[c_new(m,:,:),k_prime_new(m,:,:),y_new(m,:,:),diffs(m,:),flag(m),sol_y_new(m,:),sol_y_2_new(m,:),sol_nu_new(m,:),sol_nu_2_new(m,:),A_new(m)]...
			= OptPolUpdate_AggShock_PPA(iter,m,squeeze(grid(m,:,:,:)),...
			  k_prime,y,sol_y,Sol.sol_y_1,sol_y_2,sol_nu,sol_nu_2,A,StaticParams,Sol,Poly,...
			  squeeze(agK_prime_new_var(m,:,:,:)),squeeze(c_pr_var(m,:,:,:)));
        end		
		Sol.diff(iter,:)=[max(diffs(:,1)),max(diffs(:,2))];
		Sol.flag(iter,:)=flag;
        Sol.sol_y_1 = Sol.wealth-permute(reshape(sol_y_new,[size(Sol.wealth,1),size(Sol.wealth,3),size(Sol.wealth,2)]),[1,3,2]);
		Sol.k_prime = k_prime_new;
        Sol.c = c_new;
        Sol.y = y_new;
        Sol.sol_y = sol_y_new;
        Sol.sol_y_2 = sol_y_2_new;
        Sol.sol_nu = sol_nu_new;
        Sol.sol_nu_2 = sol_nu_2_new;
        Sol.A = A_new;
        
        iter = iter+1;
        
	end	
    Sol.time_solve = toc;
end
Sol.k_prime = Sol.k_prime.*(Sol.k_prime>eps);
if StaticParams.methodEqPFI == 0
	save(strcat(folder,'/Sol',num2str(M),'_PPA.mat'),'-struct','Sol'); 
else
	save(strcat(folder,'/Sol',num2str(M),'_PFI.mat'),'-struct','Sol'); 
end
end
end
end