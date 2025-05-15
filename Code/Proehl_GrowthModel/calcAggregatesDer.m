% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function computes the derivative w.r.t. individual capital of
% aggregate capital and weights of the polynomial chaos which are needed in
% the computation of the Jacobian in the proximal point algorithm.
%__________________________________________________________________________
function [K_der_small,w_der_small] = calcAggregatesDer(sol_k,grid,StaticParams,Poly,h)
    K_der_small = zeros(2,2,length(StaticParams.kGrid_pol));
    [~,lowIdx,upIdx] = getCapIdx(StaticParams.kGrid_pol,StaticParams.kGrid);
    gridIdx = unique([lowIdx',upIdx',(length(StaticParams.kGrid)-2:length(StaticParams.kGrid))]);
    der_grid = StaticParams.kGrid(gridIdx);
    H = zeros(length(StaticParams.kGrid),length(der_grid));
    H(gridIdx,:) = diag(h*ones(1,length(der_grid)));
    bump = zeros(2,4,length(StaticParams.kGrid),length(der_grid));
    bump(1,1,:,:) = H; %bump unemployed-bad
    bump(1,3,:,:) = H; %bump unemployed-good
    bump(2,2,:,:) = H; %bump employed-bad
    bump(2,4,:,:) = H; %bump employed-good
    K_der = zeros(2,2,length(der_grid));
    if nargout>1
		w_der_small = zeros(2,2,2,Poly.phi_N+1,length(StaticParams.kGrid_pol));
        w_der = zeros(2,2,2,Poly.phi_N+1,length(der_grid));
    end
    if size(sol_k,1)<4
       sol_k = repmat(sol_k,[2,1,1]); 
    end
    for z=1:2
        [a1,b1] = calcAggregates((repmat(sol_k,[1,1,length(der_grid)])+squeeze(bump(z,:,:,:))),grid,StaticParams,Poly);
        [a2,b2] = calcAggregates((repmat(sol_k,[1,1,length(der_grid)])-squeeze(bump(z,:,:,:))),grid,StaticParams,Poly);	
        K_der(z,:,:) = (a1-a2)/(2*h);
        K_der(z,:,end-1:end) = repmat(K_der(z,:,end-2),[1,1,2])...
                               +repmat((K_der(z,:,end-2)-K_der(z,:,end-3))./(StaticParams.kGrid_pol(end-2)-StaticParams.kGrid_pol(end-3)),[1,1,2])...
                               .*repmat(permute((StaticParams.kGrid_pol(end-1:end)-StaticParams.kGrid_pol(end-2)),[3,1,2]),[1,2,1]);
        K_der_small(z,:,:) = interp1(der_grid,squeeze(K_der(z,:,:))',StaticParams.kGrid_pol)';
        if nargout>1
            w_der(z,:,:,:,:) = (b1-b2)/(2*h);
            w_der(z,:,:,:,end-1:end) = repmat(w_der(z,:,:,:,end-2),[1,1,1,1,2])...
                           +repmat((w_der(z,:,:,:,end-2)-w_der(z,:,:,:,end-3))./(StaticParams.kGrid_pol(end-2)-StaticParams.kGrid_pol(end-3)),[1,1,1,1,2])...
                           .*repmat(permute((StaticParams.kGrid_pol(end-1:end)-StaticParams.kGrid_pol(end-2)),[3,4,5,1,2]),[1,2,2,Poly.phi_N+1,1]);
            for n=1:2
                for i=1:2
                w_der_small(z,n,i,:,:) = interp1(der_grid,squeeze(w_der(z,n,i,:,:))',StaticParams.kGrid_pol)';
                end
            end
        end
    end
end