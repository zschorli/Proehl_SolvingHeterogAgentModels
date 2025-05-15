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

% Use either the Proximal point algorithm (method = 0) or policy function
% iteration (method = 1)
for method = 1:-1:0
StaticParams.methodEqPFI = method;  

% Vary the different truncation orders of the polynomial chaos expansion.
for M=0:order
    
%--------------------------------------------------------------------------
% Set the grid for polynomial chaos coefficients
%--------------------------------------------------------------------------
if StaticParams.methodEqPFI == 0
	Poly =load(strcat(folder,'/Poly',num2str(M),'_PPA.mat'));
else
	Poly =load(strcat(folder,'/Poly',num2str(M),'_PFI.mat'));
end

Sol = struct;
switch case_nr
    case 1
        grid_N1 = 5;
        grid_N2 = 5^(Poly.phi_N>0);
        grid_N3 = 2^(Poly.phi_N>1);
        grid_N4 = 2^(Poly.phi_N>2);
        grid_N5 = 2^(Poly.phi_N>3);
        Sol.distrGrid_min = [35, 0.7, -1e-3, -3e-5, -1e-6];
        Sol.distrGrid_max = [42, 1.5,  8e-3,  1e-5,  1e-6];
        Sol.grid_N = [grid_N1, grid_N2, grid_N3, grid_N4, grid_N5];
        G = zeros(5,grid_N5,grid_N4,grid_N3,grid_N2,grid_N1);
        [G(5,:,:,:,:,:),G(4,:,:,:,:,:),G(3,:,:,:,:,:),G(2,:,:,:,:,:),G(1,:,:,:,:,:)] = ...
        ndgrid(1:Sol.grid_N(5),1:Sol.grid_N(4),1:Sol.grid_N(3),1:Sol.grid_N(2),1:Sol.grid_N(1));
    case 2
        grid_N1 = 5;
        grid_N2 = 5^(Poly.phi_N>0);
        grid_N3 = 3^(Poly.phi_N>1);
        grid_N4 = 2^(Poly.phi_N>2);
        grid_N5 = 2^(Poly.phi_N>3);
        Sol.distrGrid_min = [35, 0.7, -1e-3, -3e-5, -1e-6];
        Sol.distrGrid_max = [42, 1.5,  8e-3,  1e-5,  1e-6];
        Sol.grid_N = [grid_N1, grid_N2, grid_N3, grid_N4, grid_N5];
        G = zeros(5,grid_N5,grid_N4,grid_N3,grid_N2,grid_N1);
        [G(5,:,:,:,:,:),G(4,:,:,:,:,:),G(3,:,:,:,:,:),G(2,:,:,:,:,:),G(1,:,:,:,:,:)] = ...
        ndgrid(1:Sol.grid_N(5),1:Sol.grid_N(4),1:Sol.grid_N(3),1:Sol.grid_N(2),1:Sol.grid_N(1));
    case 3
        grid_N1 = 9;
        grid_N2 = 9^(Poly.phi_N>0);
        grid_N3 = 7^(Poly.phi_N>1);
        grid_N4 = 3^(Poly.phi_N>2);
        grid_N5 = 2^(Poly.phi_N>3);
        Sol.distrGrid_min = [35, 0.7, -1e-3, -3e-5, -1e-6];
        Sol.distrGrid_max = [42, 1.5,  8e-3,  1e-5,  1e-6];
        Sol.grid_N = [grid_N1, grid_N2, grid_N3, grid_N4, grid_N5];
        G = zeros(5,grid_N5,grid_N4,grid_N3,grid_N2,grid_N1);
        [G(5,:,:,:,:,:),G(4,:,:,:,:,:),G(3,:,:,:,:,:),G(2,:,:,:,:,:),G(1,:,:,:,:,:)] = ...
        ndgrid(1:Sol.grid_N(5),1:Sol.grid_N(4),1:Sol.grid_N(3),1:Sol.grid_N(2),1:Sol.grid_N(1));
end
grid_N = grid_N1*grid_N2*grid_N3*grid_N4*grid_N5;

Sol.distrGrid = zeros(max(Sol.grid_N),Poly.phi_N+1);
Sol.idx = zeros(grid_N,Poly.phi_N+1);
for i=1:Poly.phi_N+1
    if i>1
        Sol.distrGrid(1:Sol.grid_N(i),i) = linspace(Sol.distrGrid_min(i),Sol.distrGrid_max(i),Sol.grid_N(i))';
    else
        Sol.distrGrid(1:Sol.grid_N(i),i) = sort([linspace(Sol.distrGrid_min(i),Sol.distrGrid_max(i),Sol.grid_N(i)-2),StaticParams.kss_b,StaticParams.kss_g])';
    end
    Sol.idx(:,i) = reshape(squeeze(G(i,:,:,:,:,:)),[],1);
end
 
clear i;

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
	Sol.flag = ones(1,grid_N);
	y_new = zeros(grid_N,4,length(StaticParams.kGrid_pol));
	flag = zeros(grid_N,1);
end

grid = zeros(grid_N,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:)));
diffs = 1e8*ones(grid_N,2);

if M==0
    if StaticParams.methodEqPFI == 0
	    preSol =load(strcat(folder,'/preSol_PPA.mat'));
    else
	    preSol =load(strcat(folder,'/preSol_PFI.mat'));
    end
    c_g_approx = @(z_idx,k,K) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,squeeze(preSol.c_g(z_idx,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)));
    c_b_approx = @(z_idx,k,K) interpn(StaticParams.kGrid_pol,StaticParams.KGrid,squeeze(preSol.c_b(z_idx,:,:)),k,min(StaticParams.KGrid(end),max(StaticParams.KGrid(1),K)));
    K_ss = [(2*StaticParams.p(1:2)'*preSol.pdf_cond_b*StaticParams.kGrid');...
            (2*StaticParams.p(3:4)'*preSol.pdf_cond_g*StaticParams.kGrid')];
    c_ss = [c_b_approx(1,StaticParams.kGrid_pol,K_ss(1))';...
            c_b_approx(2,StaticParams.kGrid_pol,K_ss(1))';...
            c_g_approx(1,StaticParams.kGrid_pol,K_ss(2))';...
            c_g_approx(2,StaticParams.kGrid_pol,K_ss(2))'];
else
    if StaticParams.methodEqPFI == 0
	    preSol = load(strcat(folder,'/Sol',num2str(M-1),'_PPA.mat')); 
    else
	    preSol = load(strcat(folder,'/Sol',num2str(M-1),'_PFI.mat')); 
    end
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
    
    % interest rate 
    Sol.rate(m,:,:) = [repmat(StaticParams.rate(0,Sol.agK(m)),[1,numel(StaticParams.kGrid_pol)]);...
                       repmat(StaticParams.rate(1,Sol.agK(m)),[1,numel(StaticParams.kGrid_pol)])];
    % wage 
    Sol.wage(m,:,:) = [repmat(StaticParams.wage(0,Sol.agK(m)),[1,numel(StaticParams.kGrid_pol)]);...
                       repmat(StaticParams.wage(1,Sol.agK(m)),[1,numel(StaticParams.kGrid_pol)])];   
    % income
    Sol.wealth(m,:,:) = [StaticParams.wealth(0,Sol.agK(m),StaticParams.kGrid_pol);...
                         StaticParams.wealth(1,Sol.agK(m),StaticParams.kGrid_pol)]; 
    if M==0
        Sol.c(m,:,:) = min(squeeze(Sol.wealth(m,:,:)),c_ss);
        Sol.k_prime(m,:,:) = squeeze(Sol.wealth(m,:,:)-Sol.c(m,:,:));
    end
end
if M>0
    aux = reshape(repmat(1:size(preSol.k_prime,1),[Sol.grid_N(M+1),1]),[],1);
    Sol.k_prime = preSol.k_prime(aux,:,:);
    Sol.c = preSol.c(aux,:,:);
    if StaticParams.methodEqPFI == 0
	    Sol.y = preSol.y(aux,:,:);
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
    = OptPol_AggShock_PFI(Sol.c,Sol.k_prime,grid,StaticParams,Poly,Sol.grid_N,Sol.distrGrid,Sol.idx,Sol.wealth);
    Sol.time_solve = toc;
else
	agK_prime_new_var = zeros(grid_N,2,2,length(StaticParams.kGrid_pol));
	c_pr_var = zeros(grid_N,4,4,length(StaticParams.kGrid_pol));
    mult = max(1,M)/10;
	tic;
	while max(Sol.diff(end,:)) > StaticParams.criter_k && sum(Sol.flag(end,:)>0)>0
        parfor m=1:grid_N  %
            k_prime = squeeze(Sol.k_prime(m,:,:));
			y = squeeze(Sol.y(m,:,:));
            
            % compute the derivatives of aggregate variables w.r.t. changes
            % in individual capital (needed for the Jacobian in the
            % nonlinear solver which solves the Euler equation)
            condition = iter<=3 ...
                        || sum(sum(sum(c_pr_var(m,:,:,:),2),3),4)==0 ...
                        || floor(max(sum(Sol.diff(3:end,:),1))/norm(3*Sol.diff(1,:)))...
                           >floor(max(sum(Sol.diff(3:end-1,:),1))/norm(3*Sol.diff(1,:)));
            if condition
				[agK_prime_new_var(m,:,:,:),c_pr_var(m,:,:,:)] = calcJacobianProxies(k_prime,squeeze(grid(m,:,:,:)),StaticParams,Sol,Poly);
            end 
            
            % update of the proximal point algorithm
			[c_new(m,:,:),k_prime_new(m,:,:),y_new(m,:,:),diffs(m,:),flag(m)]...
			= OptPolUpdate_AggShock_PPA(iter,m,squeeze(grid(m,:,:,:)),...
			  k_prime,y,StaticParams,Sol,Poly,mult,...
			  squeeze(agK_prime_new_var(m,:,:,:)),squeeze(c_pr_var(m,:,:,:)));
        end		
		Sol.diff(iter,:)=[max(diffs(:,1)),max(diffs(:,2))];
		Sol.flag(iter,:)=flag;
        Sol.k_prime = k_prime_new;
        Sol.c = c_new;
        Sol.y = y_new;
        
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