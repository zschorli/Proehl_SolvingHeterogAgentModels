% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This function solves the model without aggregate risk.
%__________________________________________________________________________
function solveModel_NoAggShock(StaticParams,folder,case_nr)

% Use either the Proximal point algorithm (method = 0) or policy function
% iteration (method = 1)
for method = 1:-1:0
StaticParams.methodEqPFI = method;   

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
preSol = struct;
N = 20; 
switch case_nr
    case {5,13}
        preSol.KGrid = linspace(2,15,N); % aggregate capital grid 
    case {6,14}
        preSol.KGrid = linspace(50,70,N); % aggregate capital grid 
    case {9,17}
        preSol.KGrid = linspace(2,15,N); % aggregate capital grid 
    case {10,18}
        preSol.KGrid = linspace(85,120,N); % aggregate capital grid 
    case {11,19}
        preSol.KGrid = linspace(80,100,N); % aggregate capital grid 
    case {12,20}
        preSol.KGrid = linspace(2,15,N); % aggregate capital grid 
    otherwise
        preSol.KGrid = linspace(32,45,N); % aggregate capital grid
end

% Production function quantities___________________________________________
% interest rate 
preSol.rate = StaticParams.alpha.*...
              repmat([(1-StaticParams.delta_a)*ones(2,1);(1+StaticParams.delta_a)*ones(2,1)],[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)])...
              .*(repmat(permute(preSol.KGrid,[1,3,2]),[4,length(StaticParams.kGrid_pol),1])...
              ./(StaticParams.l_bar.*repmat([StaticParams.er(1)*ones(2,1);StaticParams.er(2)*ones(2,1)],[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)])))...
              .^(StaticParams.alpha-1);   
% wage 
preSol.wage = (1-StaticParams.alpha).*...
              repmat([(1-StaticParams.delta_a)*ones(2,1);(1+StaticParams.delta_a)*ones(2,1)],[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)])...
              .*(repmat(permute(preSol.KGrid,[1,3,2]),[4,length(StaticParams.kGrid_pol),1])...
              ./(StaticParams.l_bar.*repmat([StaticParams.er(1)*ones(2,1);StaticParams.er(2)*ones(2,1)],[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)])))...
              .^(StaticParams.alpha);   
% income
preSol.wealth = ((StaticParams.l_bar-...
                repmat([StaticParams.mu.*StaticParams.ur(1)/StaticParams.er(1)*ones(2,1);...
                        StaticParams.mu.*StaticParams.ur(2)/StaticParams.er(2)*ones(2,1)],[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)]))...
                .*repmat([0;1;0;1],[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)])...
                +StaticParams.mu.*(1-repmat([0;1;0;1],[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)])))...
                .*preSol.wage...
                +(1-StaticParams.delta+preSol.rate)...
                .*repmat(StaticParams.kGrid_pol,[4,1,length(preSol.KGrid)]);

%--------------------------------------------------------------------------
% Computing the solution when the aggregate state is fixed to 1 (good state)
%--------------------------------------------------------------------------     
c_g = repmat(log(StaticParams.kGrid_pol+1)+0.2,[2,1,N]);
k_prime_g = squeeze(preSol.wealth(3:4,:,:))-c_g;

if StaticParams.methodEqPFI == 1
	tic;
	[c_g,k_prime_g,diff_g]...
    = OptPol_NoAggShock_PFI(true,c_g,k_prime_g,StaticParams,preSol);
    preSol.time_g = toc;
else
	iter_ct = 2000;                   % maximal iteration count
    y_g = zeros(2,length(StaticParams.kGrid_pol),N);
    flag_g = zeros(iter_ct,N);
    diff_g = zeros(iter_ct,2,N);

    tic;
    parfor i=1:N %
        [c_g(:,:,i),k_prime_g(:,:,i),y_g(:,:,i),diff_g(:,:,i),flag_g(:,i)]...
        = OptPol_NoAggShock_PPA(i,true,squeeze(c_g(:,:,i)),...
          squeeze(k_prime_g(:,:,i)),squeeze(y_g(:,:,i)),StaticParams,preSol,iter_ct);
    end
    preSol.time_g = toc;
end
preSol.k_prime_g = k_prime_g.*(k_prime_g>eps);
preSol.c_g = c_g;
preSol.diffs_g = diff_g;

if StaticParams.methodEqPFI == 0
    preSol.y_g = y_g;
    preSol.flags_g = flag_g;
	save(strcat(folder,'/preSol_PPA.mat'),'-struct','preSol');
else
	save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');
end

%--------------------------------------------------------------------------
% Computing the solution when the aggregate state is fixed to 0 (bad state)
%--------------------------------------------------------------------------
c_b = repmat(log(StaticParams.kGrid_pol+1)+0.2,[2,1,N]);
k_prime_b = squeeze(preSol.wealth(1:2,:,:))-c_b;

if StaticParams.methodEqPFI == 1
	tic;
	[c_b,k_prime_b,diff_b]...
    = OptPol_NoAggShock_PFI(false,c_b,k_prime_b,StaticParams,preSol);
    preSol.time_b = toc;
else
	iter_ct = 2000;                   % maximal iteration count
	y_b = zeros(2,length(StaticParams.kGrid_pol),N);
    flag_b = zeros(iter_ct,N); 
    diff_b = zeros(iter_ct,2,N);

	tic;
	parfor i=1:N %
		[c_b(:,:,i),k_prime_b(:,:,i),y_b(:,:,i),diff_b(:,:,i),flag_b(:,i)]...
		= OptPol_NoAggShock_PPA(i,false,squeeze(c_b(:,:,i)),...
          squeeze(k_prime_b(:,:,i)),squeeze(y_b(:,:,i)),StaticParams,preSol,iter_ct);
	end
	preSol.time_b = toc;
end
preSol.k_prime_b = k_prime_b.*(k_prime_b>eps);
preSol.c_b = c_b;
preSol.diffs_b = diff_b;

if StaticParams.methodEqPFI == 0
    preSol.y_b = y_b;
	preSol.flags_b = flag_b;
	save(strcat(folder,'/preSol_PPA.mat'),'-struct','preSol');
else
	save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');
end

clear N c k_prime y k_prime_g k_prime_b c_g c_b y_g y_b diff_g diff_b ...
      flag_g flag_b;
  
%--------------------------------------------------------------------------
% Computing the stationary distributions for the cases without aggregate
% shocks
%--------------------------------------------------------------------------
pdf = zeros(2,2,length(StaticParams.kGrid));  
k = zeros([2,size(preSol.k_prime_g)]);
k(1,:,:,:) = preSol.k_prime_b;
k(2,:,:,:) = preSol.k_prime_g;
    
tic;
parfor bool_gb=0:1 %
    k_prime = max(0,squeeze(k(bool_gb+1,:,:,:)));
    pdf_cond = zeros(2,length(StaticParams.kGrid));  
    pdf_cond(:,getCapIdx(StaticParams.kss,StaticParams.kGrid)) = 1;
    k_prime_approx = @(z_idx,k,K) interpn(StaticParams.kGrid_pol,preSol.KGrid,squeeze(k_prime(z_idx,:,:)),k,K);
    iter = 1;
    diffs = 1e8*ones(2,1);
    while (diffs(end)) > StaticParams.criter_pdf
        pdf_cond_new = (repmat(1./(StaticParams.p(2*bool_gb+(1:2))'*StaticParams.P(2*bool_gb+(1:2),2*bool_gb+(1:2)))',[1,2])...
                       .*StaticParams.P(2*bool_gb+(1:2),2*bool_gb+(1:2))')...
                       *(repmat(StaticParams.p(2*bool_gb+(1:2)),[1,size(pdf_cond,2)])...
                       .*pdf_cond);
        cdf_cond_new = cummax(min(1,cumsum(pdf_cond_new,2)),2);
        cdf_cond_new(1,cdf_cond_new(1,:)==max(cdf_cond_new(1,:))) = 1;
        cdf_cond_new(2,cdf_cond_new(2,:)==max(cdf_cond_new(2,:))) = 1;
        agK = StaticParams.kGrid*pdf_cond_new'*[StaticParams.ur(bool_gb+1);StaticParams.er(bool_gb+1)];
        k_prime = sort([k_prime_approx(1,StaticParams.kGrid_pol,min(preSol.KGrid(end),max(preSol.KGrid(1),agK)))';...
                        k_prime_approx(2,StaticParams.kGrid_pol,min(preSol.KGrid(end),max(preSol.KGrid(1),agK)))'],2);
        pdf_cond_new = calcEoPDistr(max(StaticParams.k_min,interp1(StaticParams.kGrid_pol,k_prime',StaticParams.kGrid)'),...
                       cdf_cond_new,StaticParams.kGrid);
                
        diffs(iter) = norm(pdf_cond_new-pdf_cond,2);
        pdf_cond_new = pdf_cond_new.*(pdf_cond_new>eps);
        pdf_cond = pdf_cond_new./repmat(sum(pdf_cond_new,2),[1,size(pdf_cond_new,2)]);
        iter = iter+1;
    end
    pdf(bool_gb+1,:,:) = pdf_cond;
end
preSol.time_distr = toc;

preSol.pdf_cond_b = squeeze(pdf(1,:,:)); 
preSol.pdf_cond_g = squeeze(pdf(2,:,:)); 
preSol.cdf_cond_b = cummax(min(1,cumsum(preSol.pdf_cond_b,2)),2); 
preSol.cdf_cond_g = cummax(min(1,cumsum(preSol.pdf_cond_g,2)),2);
preSol.cdf_cond_b(:,end) = 1;
preSol.cdf_cond_g(:,end) = 1;
clear k k_prime k_prime_approx pdf pdf_cond pdf_cond_new cdf_cond agK ...
      diffs iter bool_gb i iter_ct;
  
if StaticParams.methodEqPFI == 0
	save(strcat(folder,'/preSol_PPA.mat'),'-struct','preSol');
else
	save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');
end
end
end