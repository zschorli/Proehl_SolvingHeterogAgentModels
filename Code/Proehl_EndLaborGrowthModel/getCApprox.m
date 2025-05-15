% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function computes the approximated consumption (interpolation w.r.t.
% the polynomial chaos weights) 
%__________________________________________________________________________
function c = getCApprox(grid_N,grid,M_red,M_red_idx,x)

M = zeros([size(M_red,2),grid_N(1:size(grid,2))]);
for m=1:size(M_red_idx,1)
    d = num2cell(M_red_idx(m,:)); 
    M(:,d{:}) = M_red(m,:);
end

a = num2cell(repmat((1:2)',[1,size(grid,2)]),1);
A = zeros([size(grid,2),2*ones(1,size(grid,2))]);
b=[];
for i=1:size(grid,2)
    idx = max(1,min(size(grid,1)-1-sum((grid(:,i)==0)),sum(x(i)>grid(logical(1-(grid(:,i)==0)),i))));
    b(:,end+1)=[idx;idx+1];
    A(i,a{:}) = repmat(permute([grid(idx+1,i)-x(i);x(i)-grid(idx,i)]./(grid(idx+1,i)-grid(idx,i)),[(3:2+i-1),1:2]),[2*ones(1,i-1),1,2*ones(1,length(i+1:size(grid,2)))]);
end
b = num2cell(b,1);
A = repmat(prod(A,1),[size(M,1),ones(1,size(grid,2))]);
c = A.*M(:,b{:});
for i=1:size(grid,2)
c = sum(c,1+i);
end
c = squeeze(c)';

end

