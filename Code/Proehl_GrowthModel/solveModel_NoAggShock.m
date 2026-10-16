% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
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
n = length(StaticParams.kGrid_pol); 
N = length(StaticParams.KGrid); 

%--------------------------------------------------------------------------
% Computing the solution when the aggregate state is fixed to 1 (good state)
%--------------------------------------------------------------------------     
c_g = repmat(log(StaticParams.kGrid_pol+1)+0.2,[2,1,N]);
k_prime_g =(StaticParams.wealth(1,StaticParams.KGrid,StaticParams.kGrid_pol)-c_g);

if StaticParams.methodEqPFI == 1
	tstart=cputime;
	[c_g,k_prime_g,diff_g]...
    = OptPol_NoAggShock_PFI(true,c_g,k_prime_g,StaticParams);
    preSol.time_g = cputime-tstart;
else
	iter_ct = 5000;                   % maximal iteration count
    y_g = zeros(2,n,N);
    flag_g = zeros(iter_ct,N);
    diff_g = zeros(iter_ct,2,N);

    tstart=cputime;
    parfor i=1:N %
        [c_g(:,:,i),k_prime_g(:,:,i),y_g(:,:,i),diff_g(:,:,i),flag_g(:,i)]...
        = OptPol_NoAggShock_PPA(i,1,squeeze(k_prime_g(:,:,i)),...
          squeeze(y_g(:,:,i)),StaticParams,iter_ct);
    end
    preSol.time_g = cputime-tstart;
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
k_prime_b = StaticParams.wealth(0,StaticParams.KGrid,StaticParams.kGrid_pol)-c_b;

if StaticParams.methodEqPFI == 1
	tstart=cputime;
	[c_b,k_prime_b,diff_b]...
    = OptPol_NoAggShock_PFI(false,c_b,k_prime_b,StaticParams);
    preSol.time_b = cputime-tstart;
else
	iter_ct = 5000;                   % maximal iteration count
	y_b = zeros(2,n,N);
    flag_b = zeros(iter_ct,N); 
    diff_b = zeros(iter_ct,2,N);

	tstart=cputime;
	parfor i=1:N %
		[c_b(:,:,i),k_prime_b(:,:,i),y_b(:,:,i),diff_b(:,:,i),flag_b(:,i)]...
		= OptPol_NoAggShock_PPA(i,0,squeeze(k_prime_b(:,:,i)),...
          squeeze(y_b(:,:,i)),StaticParams,iter_ct);
	end
	preSol.time_b = cputime-tstart;
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
P = zeros([2,size(StaticParams.P_g)]);
P(1,:,:) = StaticParams.P_b;
P(2,:,:) = StaticParams.P_g;   

tstart=cputime;
parfor bool_gb=0:1 %
    k_prime = max(0,squeeze(k(bool_gb+1,:,:,:)));

    [Yz,Yk] = ndgrid(1:2,StaticParams.kGrid_pol,1);
    k_prime_approx = @(K) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,k_prime,Yz,Yk,K.*ones(size(Yz)));
    p = StaticParams.p(2*bool_gb+(1:2))*2;
    agK_stationary = fzero(@(K) K-calcStationaryDistr(k_prime_approx(K),squeeze(P(bool_gb+1,:,:)),p,StaticParams),[35,42]);
    [~,pdf_stationary] = calcStationaryDistr(k_prime_approx(agK_stationary),squeeze(P(bool_gb+1,:,:)),p,StaticParams);
    
    pdf(bool_gb+1,:,:) = pdf_stationary;
end
preSol.time_distr = cputime-tstart;

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