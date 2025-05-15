% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file computes the stationary distribution for given policy
% functions
%__________________________________________________________________________
function [agK,agL,pdf,pdf_multistep] = calcStationaryDistr(k_prime,l,P,p,StaticParams,pdf_old,steps)
% pdf is function of (z^id,k)
N = length(p);
kGrid = StaticParams.kGrid;
pdf = zeros(N,length(kGrid));
pdf_large = reshape(pdf,[],1);

%--------------------------------------------------------------------------
% Find the surrounding grid points of the new policy values
%--------------------------------------------------------------------------
[Yz,Yk] = ndgrid(1:N,StaticParams.kGrid);
k_int = interpn(1:N,StaticParams.kGrid_pol,k_prime,Yz,Yk);
l_int = interpn(1:N,StaticParams.kGrid_pol,l,Yz,Yk);

% last argument: 1=lower grid idx, 2=upper grid idx     
k_bound = zeros([size(Yk),2]);
for i=1:length(kGrid)
    k_bound(:,:,1) = k_bound(:,:,1) + (k_int>kGrid(i));
end
aux = squeeze(k_bound(:,:,1));
aux(k_int>StaticParams.k_min) = min(aux(k_int>StaticParams.k_min)+1,length(kGrid));
k_bound(:,:,2) = aux;
k_bound(k_bound==0) = 1; 
k_idx = zeros(size(k_bound));
for i=1:N
k_idx(i,:,:) = sub2ind([N,numel(kGrid)],i.*ones(size(k_bound(i,:,:))),k_bound(i,:,:));
end

%--------------------------------------------------------------------------
% Find the weights for the surrounding grid values (by proximity)
%--------------------------------------------------------------------------
% histogram convention: probability assigned to a grid point is the
% probability that savings are <= the grid point but > the grid point to 
% the left
k_weight = zeros(size(k_bound));
low = squeeze(k_bound(:,:,1));
up = squeeze(k_bound(:,:,2));
idx = sub2ind([N,length(kGrid)],repmat((1:N)',[1,length(kGrid)]),repmat((1:length(kGrid)),[N,1]));
idx_1 = sub2ind([N,length(kGrid)],repmat((1:N)',[1,length(kGrid)]),repmat([1,1:length(kGrid)-1],[N,1]));
weight_up = zeros(size(up));
weight_up(low<up) = (k_int(idx(low<up))-kGrid(low(low<up))')...
                  ./(k_int(idx(low<up))-min(kGrid(low(low<up))',k_int(idx_1(low<up))));
weight_up(low==up) = 0.5;
k_weight(:,:,2) = weight_up;
k_weight(:,:,1) = 1-weight_up;

clear low up weight_up;

%--------------------------------------------------------------------------
% Construct the transition matrix for the pdf of current holdings to
% next-period holdings
%--------------------------------------------------------------------------
% first compute the transition from EoP distribution to BoP distribution
% due to the realization of the exogeneous shocks
P_aux = repmat(1./(p'*P)',[1,N]).*P'.*repmat(p',[N,1]);
idx_a = zeros(N,N,numel(kGrid));
idx_b = zeros(N,N,numel(kGrid));
for i=1:N
    idx_a(i,:,:) = repmat(sub2ind([N,numel(kGrid)],i.*ones(1,numel(kGrid)),1:numel(kGrid)),[N,1]);
    idx_b(:,i,:) = repmat(sub2ind([N,numel(kGrid)],i.*ones(1,numel(kGrid)),1:numel(kGrid)),[N,1]);
end
val = repmat(P_aux,[1,1,numel(kGrid)]);
pdf_trans1 = sparse(reshape(idx_a,[],1),reshape(idx_b,[],1),reshape(val,[],1),numel(pdf_large),numel(pdf_large));

% then compute the transition due to the optimal choices of the agents
pdf_idx = (1:numel(pdf_large))';
val = zeros(numel(pdf_large)^2,1);
for k=1:2
    pdf_idx_pr = reshape(squeeze(k_idx(:,:,k)),[],1);
    weight = reshape(squeeze(k_weight(:,:,k)),[],1);
    idx = sub2ind([numel(pdf_large),numel(pdf_large)],pdf_idx_pr,pdf_idx);
    val(idx) = val(idx)+weight;
end
[idx_a,idx_b] = ind2sub([numel(pdf_large),numel(pdf_large)],1:numel(pdf_large)^2);
pdf_trans2 = sparse(idx_a,idx_b,val,numel(pdf_large),numel(pdf_large));

pdf_trans = pdf_trans2*pdf_trans1;

%--------------------------------------------------------------------------
% Compute the fixed point of the transition and aggregate holdings
%--------------------------------------------------------------------------
[V,~] = eigs(pdf_trans,1,1);
if abs(min(V))>max(V)
    V=-V;
end
pdf_large = V;
pdf_large(pdf_large<=eps) = 0;
pdf = reshape(pdf_large,size(pdf));
pdf = pdf./sum(pdf,2);
agK = sum(pdf.*Yk,2)'*p;
p_aux = zeros(size(p));
p_aux(2:2:numel(p_aux)) = p(2:2:numel(p_aux));
p_aux = p_aux./sum(p_aux);
agL = sum(pdf.*l_int,2)'*p_aux;
if nargin>=6
    pdf_large = reshape(pdf_old,[],1);
    for i=1:steps
    pdf_large = pdf_trans*pdf_large;
    end
    pdf_large(pdf_large<=eps) = 0;
    pdf_multistep = reshape(pdf_large,size(pdf));
    pdf_multistep = pdf_multistep./sum(pdf_multistep,2);
end
end