% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This file computes the stationary distribution for given policy
% functions
%__________________________________________________________________________
function [agA,agH,pdf,pdf_multistep] = calcStationaryDistr(a_prime,h,P,p,StaticParams,pdf_old,steps)
% pdf is function of (z^id,a)
N = length(p);
aGrid = StaticParams.aGrid;
pdf = zeros(N,length(aGrid));
pdf_large = reshape(pdf,[],1);

%--------------------------------------------------------------------------
% Find the surrounding grid points of the new policy values
%--------------------------------------------------------------------------
[Yz,Ya] = ndgrid(1:N,StaticParams.aGrid);
a_int = interpn(1:N,StaticParams.aGrid_pol,a_prime,Yz,Ya);
h_int = interpn(1:N,StaticParams.aGrid_pol,h,Yz,Ya);

% last argument: 1=lower grid idx, 2=upper grid idx     
a_bound = zeros([size(Ya),2]);
for i=1:length(aGrid)
    a_bound(:,:,1) = a_bound(:,:,1) + (a_int>aGrid(i));
end
aux = squeeze(a_bound(:,:,1));
aux(a_int>StaticParams.a_min) = min(aux(a_int>StaticParams.a_min)+1,length(aGrid));
a_bound(:,:,2) = aux;
a_bound(a_bound==0) = 1; 
a_idx = zeros(size(a_bound));
for i=1:N
a_idx(i,:,:) = sub2ind([N,numel(aGrid)],i.*ones(size(a_bound(i,:,:))),a_bound(i,:,:));
end

%--------------------------------------------------------------------------
% Find the weights for the surrounding grid values (by proximity)
%--------------------------------------------------------------------------
% histogram convention: probability assigned to a grid point is the
% probability that savings are <= the grid point but > the grid point to 
% the left
a_weight = zeros(size(a_bound));
low = squeeze(a_bound(:,:,1));
up = squeeze(a_bound(:,:,2));
idx = sub2ind([N,length(aGrid)],repmat((1:N)',[1,length(aGrid)]),repmat((1:length(aGrid)),[N,1]));
idx_1 = sub2ind([N,length(aGrid)],repmat((1:N)',[1,length(aGrid)]),repmat([1,1:length(aGrid)-1],[N,1]));
weight_up = zeros(size(up));
weight_up(low<up) = (a_int(idx(low<up))-aGrid(low(low<up))')...
                  ./(a_int(idx(low<up))-min(aGrid(low(low<up))',a_int(idx_1(low<up))));
weight_up(low==up) = 0.5;
a_weight(:,:,2) = weight_up;
a_weight(:,:,1) = 1-weight_up;

clear low up weight_up;

%--------------------------------------------------------------------------
% Construct the transition matrix for the pdf of current holdings to
% next-period holdings
%--------------------------------------------------------------------------
% first compute the transition from EoP distribution to BoP distribution
% due to the realization of the exogeneous shocks
P_aux = repmat(1./(p'*P)',[1,N]).*P'.*repmat(p',[N,1]);
idx_a = zeros(N,N,numel(aGrid));
idx_b = zeros(N,N,numel(aGrid));
for i=1:N
    idx_a(i,:,:) = repmat(sub2ind([N,numel(aGrid)],i.*ones(1,numel(aGrid)),1:numel(aGrid)),[N,1]);
    idx_b(:,i,:) = repmat(sub2ind([N,numel(aGrid)],i.*ones(1,numel(aGrid)),1:numel(aGrid)),[N,1]);
end
val = repmat(P_aux,[1,1,numel(aGrid)]);
pdf_trans1 = sparse(reshape(idx_a,[],1),reshape(idx_b,[],1),reshape(val,[],1),numel(pdf_large),numel(pdf_large));

% then compute the transition due to the optimal choices of the agents
pdf_idx = (1:numel(pdf_large))';
val = zeros(numel(pdf_large)^2,1);
for a=1:2
    pdf_idx_pr = reshape(squeeze(a_idx(:,:,a)),[],1);
    weight = reshape(squeeze(a_weight(:,:,a)),[],1);
    idx = sub2ind([numel(pdf_large),numel(pdf_large)],pdf_idx_pr,pdf_idx);
    val(idx) = val(idx)+weight;
end
[idx_a,idx_b] = ind2sub([numel(pdf_large),numel(pdf_large)],1:numel(pdf_large)^2);
pdf_trans2 = sparse(idx_a,idx_b,val,numel(pdf_large),numel(pdf_large));

pdf_trans = pdf_trans2*pdf_trans1;
% pdf_trans = balance(full(pdf_trans),'noperm');

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
agA = sum(pdf.*Ya,2)'*p;
agH = sum(pdf.*h_int,2)'*p;
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