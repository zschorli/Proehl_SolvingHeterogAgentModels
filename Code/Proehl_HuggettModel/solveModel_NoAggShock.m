% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE October 2018
%
% DESCRIPTION
% This function solves the model without aggregate risk.
%__________________________________________________________________________
function solveModel_NoAggShock(StaticParams,folder)

% Use either the Proximal point algorithm (method = 0) or policy function
% iteration (method = 1)
for method = 1:-1:0
StaticParams.methodEqPFI = method;   

%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------
preSol = struct;
N = 20; 
preSol.KGrid = linspace(0.85,1.2,N); % aggregate capital grid

% Income quantities________________________________________________________
% endowment
preSol.wage = repmat([0.982*[2/3.06;4.12/3.06];...
                      1.054*[2/3.06;4.12/3.06]],[1,length(StaticParams.kGrid_pol),length(preSol.KGrid)]);
% income
preSol.wealth = preSol.wage...
                +repmat(StaticParams.kGrid_pol,[4,1,length(preSol.KGrid)]);
            
%--------------------------------------------------------------------------
% Computing the solution when the aggregate state is fixed to 1 (good state)
%--------------------------------------------------------------------------     
c_g = squeeze(preSol.wage(3:4,:,:))...
      +0.4*repmat(StaticParams.kGrid_pol,[2,1,length(preSol.KGrid)]);
k_prime_g = max(StaticParams.k_min,(squeeze(preSol.wealth(3:4,:,:))-c_g)./repmat(permute(preSol.KGrid,[3,1,2]),[2,length(StaticParams.kGrid_pol),1]));

if StaticParams.methodEqPFI == 1
	tic;
	[c_g,k_prime_g,diff_g]...
    = OptPol_NoAggShock_PFI(true,c_g,k_prime_g,StaticParams,preSol);
    preSol.time_g = toc;
else
	iter_ct = 2000;                   % maximal iteration count
    y_g = zeros(2,length(StaticParams.kGrid_pol),N);
    y2_g = zeros(2,length(StaticParams.kGrid_pol),N);
    flag_g = zeros(iter_ct,N);
    diff_g = zeros(iter_ct,3,N);

    tic;
    parfor i=1:N %
        [c_g(:,:,i),k_prime_g(:,:,i),y_g(:,:,i),y2_g(:,:,i),diff_g(:,:,i),flag_g(:,i)]...
        = OptPol_NoAggShock_PPA(i,true,squeeze(c_g(:,:,i)),...
          squeeze(k_prime_g(:,:,i)),squeeze(y_g(:,:,i)),squeeze(y2_g(:,:,i)),StaticParams,preSol,iter_ct);
    end
    preSol.time_g = toc;
end
preSol.k_prime_g = k_prime_g;
preSol.c_g = c_g;
preSol.diffs_g = diff_g;

if StaticParams.methodEqPFI == 0
    preSol.y_g = y_g;
    preSol.y2_g = y2_g;
    preSol.flags_g = flag_g;
	save(strcat(folder,'/preSol_PPA.mat'),'-struct','preSol');
else
	save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');
end

%--------------------------------------------------------------------------
% Computing the solution when the aggregate state is fixed to 0 (bad state)
%--------------------------------------------------------------------------
c_b = squeeze(preSol.wage(1:2,:,:))...
      +0.2*repmat(StaticParams.kGrid_pol,[2,1,length(preSol.KGrid)]);
k_prime_b = (squeeze(preSol.wealth(1:2,:,:))-c_b)./repmat(permute(preSol.KGrid,[3,1,2]),[2,length(StaticParams.kGrid_pol),1]);

if StaticParams.methodEqPFI == 1
	tic;
	[c_b,k_prime_b,diff_b]...
    = OptPol_NoAggShock_PFI(false,c_b,k_prime_b,StaticParams,preSol);
    preSol.time_b = toc;
else
	iter_ct = 2000;                   % maximal iteration count
	y_b = zeros(2,length(StaticParams.kGrid_pol),N);
	y2_b = zeros(2,length(StaticParams.kGrid_pol),N);
    flag_b = zeros(iter_ct,N); 
    diff_b = zeros(iter_ct,3,N);

	tic;
	parfor i=1:N %
		[c_b(:,:,i),k_prime_b(:,:,i),y_b(:,:,i),y2_b(:,:,i),diff_b(:,:,i),flag_b(:,i)]...
		= OptPol_NoAggShock_PPA(i,false,squeeze(c_b(:,:,i)),...
          squeeze(k_prime_b(:,:,i)),squeeze(y_b(:,:,i)),squeeze(y2_b(:,:,i)),StaticParams,preSol,iter_ct);
	end
	preSol.time_b = toc;
end
preSol.k_prime_b = k_prime_b;
preSol.c_b = c_b;
preSol.diffs_b = diff_b;

if StaticParams.methodEqPFI == 0
    preSol.y_b = y_b;
    preSol.y2_b = y2_b;
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
p = zeros(2,1);
    
tic;
parfor bool_gb=0:1 
    k_prime = max(StaticParams.k_min,squeeze(k(bool_gb+1,:,:,:)));
    pdf_cond = zeros(2,length(StaticParams.kGrid));  
    pdf_cond(:,getCapIdx(0,StaticParams.kGrid)) = 1;
    p(bool_gb+1) = fzero(@(p) AggrAssetAlloc(p,k_prime,pdf_cond,StaticParams,preSol.KGrid,bool_gb,true),preSol.KGrid([1,end]));
    [~,pdf_cond] = AggrAssetAlloc(p(bool_gb+1),k_prime,pdf_cond,StaticParams,preSol.KGrid,bool_gb,true);
    pdf(bool_gb+1,:,:) = pdf_cond;
end
preSol.time_distr = toc;

preSol.price_b = p(1);
preSol.price_g = p(2);
preSol.pdf_cond_b = squeeze(pdf(1,:,:)); 
preSol.pdf_cond_g = squeeze(pdf(2,:,:)); 
preSol.cdf_cond_b = cummax(min(1,cumsum(preSol.pdf_cond_b,2)),2); 
preSol.cdf_cond_g = cummax(min(1,cumsum(preSol.pdf_cond_g,2)),2);
preSol.cdf_cond_b(:,end) = 1;
preSol.cdf_cond_g(:,end) = 1;
clear k k_prime k_prime_approx pdf pdf_cond pdf_cond_new cdf_cond agK diffs iter bool_gb i iter_ct;
if StaticParams.methodEqPFI == 0
	save(strcat(folder,'/preSol_PPA.mat'),'-struct','preSol');
else
	save(strcat(folder,'/preSol_PFI.mat'),'-struct','preSol');
end
end
end