% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This function solves the model with aggregate risk.
%__________________________________________________________________________
function solveModel_AggShock(StaticParams,case_nr,folder,order)

% Use either the Proximal point algorithm (method = 0) or policy function
% iteration (method = 1)
for method =1:-1:0
StaticParams.methodEqPFI = method;  

% Vary the different truncation orders of the polynomial chaos expansion.
% Start with order 1 as order zero is fixed due to the zero net supply
% condition.
for M=1:order
    
%--------------------------------------------------------------------------
% Set the grid for polynomial chaos coefficients
%--------------------------------------------------------------------------
if StaticParams.methodEqPFI == 0
	Poly =load(strcat(folder,'/Poly',num2str(M),'_PPA.mat'));
else
	Poly =load(strcat(folder,'/Poly',num2str(M),'_PFI.mat'));
end
grid_N1 = 1;
grid_N2 = 5^(Poly.phi_N>0);
grid_N3 = 4^(Poly.phi_N>1);
grid_N4 = 3^(Poly.phi_N>2); 
grid_N = grid_N1*grid_N2*grid_N3*grid_N4^(Poly.phi_N-2);

Sol = struct;
switch case_nr
    case {1,4,5,7,9}
        Sol.distrGrid_min = [0, 0.1, -0.02, -5e-3, -5e-3];
        Sol.distrGrid_max = [0, 0.57, 0.04,  0.05,  0.02];
    case {2,10,11,12,13}
        Sol.distrGrid_min = [0, 0.01, 0,    -0.08, -0.06];
        Sol.distrGrid_max = [0, 0.45, 0.15,  0.3,   0.35];
    case 3
        Sol.distrGrid_min = [0, 0.1, 0,    -0.02, -0.02];
        Sol.distrGrid_max = [0, 0.6, 0.06,  0.15,  0.05];
    case {6,8}
        Sol.distrGrid_min = [0, 0.2,  -0.02, -0.01, -5e-3];
        Sol.distrGrid_max = [0, 0.61,  0.04,  0.04,  0.02];
end
Sol.distrGrid = zeros(max([2,grid_N1,grid_N2,grid_N3,grid_N4]),Poly.phi_N+1);
Sol.distrGrid(1:max(2,grid_N1),1) = linspace(Sol.distrGrid_min(1),Sol.distrGrid_max(1),max(2,grid_N1))';  
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
Sol.p = zeros(grid_N,2);
Sol.wage = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.wealth = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.c = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.k_prime = zeros(grid_N,4,length(StaticParams.kGrid_pol));
Sol.diff = 1e8*ones(1,3);

c_new = zeros(grid_N,4,length(StaticParams.kGrid_pol));
k_prime_new = zeros(grid_N,4,length(StaticParams.kGrid_pol));
p_new = zeros(grid_N,2);

if StaticParams.methodEqPFI == 0
	Sol.y = zeros(grid_N,4,length(StaticParams.kGrid_pol));
    Sol.y2 = zeros(grid_N,4,length(StaticParams.kGrid_pol));
	Sol.sol_y = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	Sol.sol_nu = zeros(grid_N,4*length(StaticParams.kGrid_pol));
    Sol.sol_y_1 = zeros(grid_N,4,length(StaticParams.kGrid_pol));
	Sol.sol_y_2 = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	Sol.sol_nu_2 = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	Sol.sol_y_3 = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	Sol.sol_nu_3 = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	Sol.A = 10*ones(1,grid_N);
	Sol.flag = ones(1,grid_N);

	y_new = zeros(grid_N,4,length(StaticParams.kGrid_pol));
    y_new2 = zeros(grid_N,4,length(StaticParams.kGrid_pol));
	sol_y_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	sol_nu_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	sol_y_2_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
    sol_y_3_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	sol_nu_2_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
    sol_nu_3_new = zeros(grid_N,4*length(StaticParams.kGrid_pol));
	A_new = Sol.A;
	flag = zeros(grid_N,1);
end

grid = zeros(grid_N,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:)),length(Poly.xi{3}(1,:)));
diffs = 1e8*ones(grid_N,3);

if StaticParams.methodEqPFI == 0
	preSol =load(strcat(folder,'/preSol_PPA.mat'));
else
	preSol =load(strcat(folder,'/preSol_PFI.mat'));
end
c_g_approx = @(z_idx,k,K) interpn(StaticParams.kGrid_pol,preSol.KGrid,squeeze(preSol.c_g(z_idx,:,:)),k,min(preSol.KGrid(end),max(preSol.KGrid(1),K)));
c_b_approx = @(z_idx,k,K) interpn(StaticParams.kGrid_pol,preSol.KGrid,squeeze(preSol.c_b(z_idx,:,:)),k,min(preSol.KGrid(end),max(preSol.KGrid(1),K)));
c_ss = sort([[c_b_approx(1,StaticParams.kGrid_pol,preSol.price_b)';...
             c_b_approx(2,StaticParams.kGrid_pol,preSol.price_b)'];...
             [c_g_approx(1,StaticParams.kGrid_pol,preSol.price_g)';...
             c_g_approx(2,StaticParams.kGrid_pol,preSol.price_g)']],2);  
for m=1:grid_N
    d = num2cell(Sol.idx(m,:)); 
    weights = diag(Sol.distrGrid([d{:}],:))';
    grid(m,:,:,:,:) = squeeze(sum(repmat(weights',[1,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:))]).*Poly.total,1));
           
    Sol.agK(m) = sum(sum(sum(sum(Poly.total_pdf.*squeeze(grid(m,:,:,:,:)),4),3),2),1);    
    Sol.p(m,:) = [preSol.price_b,preSol.price_g];
       
    % wage 
    Sol.wage(m,:,:) = repmat([0.982*[2/3.06;4.12/3.06];...
                      1.054*[2/3.06;4.12/3.06]],[1,length(StaticParams.kGrid_pol)]);
    % income
    Sol.wealth(m,:,:) = squeeze(Sol.wage(m,:,:))...
                +repmat(StaticParams.kGrid_pol,[4,1]);
    Sol.c(m,:,:) = c_ss;
    Sol.k_prime(m,:,:) = squeeze(Sol.wealth(m,:,:)-Sol.c(m,:,:))./repmat(Sol.p(m,[1 1 2 2])',[1,length(StaticParams.kGrid_pol)]);
    Sol.c(m,:,:) = Sol.wealth(m,:,:)- Sol.k_prime(m,:,:).*repmat(Sol.p(m,[1 1 2 2]),[1,1,length(StaticParams.kGrid_pol)]);
    if StaticParams.methodEqPFI == 0
		Sol.y(m,:,:,:) = 0;
        Sol.y2(m,:,:,:) = 0;
		Sol.sol_y(m,:,:) = reshape(squeeze(Sol.k_prime(m,:,:))',[4*length(StaticParams.kGrid_pol),1]);
		Sol.sol_nu(m,:,:) = Sol.sol_y(m,:,:);
		Sol.sol_y_1(m,:,:) = Sol.c(m,:,:);
		Sol.sol_y_2(m,:,:) = 0;
		Sol.sol_nu_2(m,:,:) = Sol.sol_y_2(m,:,:);
		Sol.sol_y_3(m,:,:) = 0;
		Sol.sol_nu_3(m,:,:) = Sol.sol_y_3(m,:,:);
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
	[Sol.c,Sol.k_prime,Sol.p,Sol.diff]...
    = OptPol_AggShock_PFI(Sol.p,Sol.c,Sol.k_prime,grid,StaticParams,Poly,Sol.distrGrid,Sol.idx,Sol.wealth);
    Sol.time_solve = toc;
else
	c_pr_var = zeros(grid_N,4,4,length(StaticParams.kGrid_pol));
    
	tic;
    while (max(Sol.diff(end,:)) > StaticParams.criter_k && sum(Sol.flag(end,:)>0)>0) || sum(Sol.flag(end,:)-10)==0
        parfor m=1:grid_N %
            k_prime = squeeze(Sol.k_prime(m,:,:));
            y = squeeze(Sol.y(m,:,:,:));
            y2 = squeeze(Sol.y2(m,:,:,:));
            sol_y = squeeze(Sol.sol_y(m,:,:));
            sol_y_2 = squeeze(Sol.sol_y_2(m,:,:));
            sol_y_3 = squeeze(Sol.sol_y_3(m,:,:));
            sol_nu = squeeze(Sol.sol_nu(m,:,:));
            sol_nu_2 = squeeze(Sol.sol_nu_2(m,:,:));
            sol_nu_3 = squeeze(Sol.sol_nu_3(m,:,:));
            A = Sol.A(m);
            pr = Sol.p(m,:);

            % compute the derivatives of aggregate variables w.r.t. changes
            % in individual capital (needed for the Jacobian in the
            % nonlinear solver which solves the Euler equation)
            condition = iter<=3 ...
                        || sum(sum(sum(c_pr_var(m,:,:,:),2),3),4)==0 ...
                        || floor(max(sum(Sol.diff(3:end,:),1))/norm(0.5*Sol.diff(1,:)))...
                           >floor(max(sum(Sol.diff(3:end-1,:),1))/norm(0.5*Sol.diff(1,:)));
            if condition
                [~,c_pr_var(m,:,:,:)] = calcJacobianProxies(k_prime,squeeze(grid(m,:,:,:,:)),StaticParams,Sol,Poly);
            end 

            % update of the proximal point algorithm
            [c_new(m,:,:,:),k_prime_new(m,:,:,:),y_new(m,:,:,:),y_new2(m,:,:,:),p_new(m,:),diffs(m,:),flag(m),sol_y_new(m,:,:),sol_y_2_new(m,:,:),sol_y_3_new(m,:,:),sol_nu_new(m,:,:),sol_nu_2_new(m,:,:),sol_nu_3_new(m,:,:),A_new(m)]...
            = OptPolUpdate_AggShock_PPA(pr,iter,m,squeeze(grid(m,:,:,:,:)),...
              y,y2,sol_y,Sol.sol_y_1,sol_y_2,sol_y_3,sol_nu,sol_nu_2,sol_nu_3,A,StaticParams,Sol,Poly,...
              squeeze(c_pr_var(m,:,:,:))); 
        end		
        Sol.diff(iter,:)=[max(diffs(:,1)),max(diffs(:,2)),max(diffs(:,3))];
        Sol.flag(iter,:)=flag;
        Sol.k_prime = k_prime_new;
        Sol.c = c_new;
        Sol.y = y_new;
        Sol.y2 = y_new2;
        Sol.p = p_new;
        Sol.sol_y = sol_y_new;
        Sol.sol_y_1 = Sol.c;
        Sol.sol_y_2 = sol_y_2_new;
        Sol.sol_y_3 = sol_y_3_new;
        Sol.sol_nu = sol_nu_new;
        Sol.sol_nu_2 = sol_nu_2_new;
        Sol.sol_nu_3 = sol_nu_3_new;
        Sol.A = A_new;

        iter = iter+1;
    end	
    Sol.time_solve = toc;
end
if StaticParams.methodEqPFI == 0
	save(strcat(folder,'/Sol',num2str(M),'_PPA.mat'),'-struct','Sol'); 
else
	save(strcat(folder,'/Sol',num2str(M),'_PFI.mat'),'-struct','Sol'); 
end
end
end
end