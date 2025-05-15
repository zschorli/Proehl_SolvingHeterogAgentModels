% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This function solves the model without aggregate risk.
%__________________________________________________________________________

folder = strcat('compareDistr_NoAggShock');        
if exist(folder,'dir') ~= 7
    mkdir(folder);
end
StaticParams = struct;
StaticParams = setModelParams(StaticParams,1);

compareSol = struct;
compareSol.muN = 22;
muGrid = [linspace(0,0.95,compareSol.muN-2),0.965,0.98];
pdf_k_g = zeros(length(muGrid),length(StaticParams.kGrid));
pdf_k_b = zeros(length(muGrid),length(StaticParams.kGrid));
pdf_c_g = zeros(length(muGrid),1001);
pdf_c_b = zeros(length(muGrid),1001);
quantileGrid = [0.01,0.05,0.1:0.1:0.9,0.95,0.99];
quantiles_b = zeros(length(muGrid),length(quantileGrid));
quantiles_g = zeros(length(muGrid),length(quantileGrid));
cdf_c_grid = linspace(0,10,1001);
fullInsurance_K_g = zeros(length(muGrid)+1,1);
fullInsurance_C_g = zeros(length(muGrid)+1,1);
fullInsurance_K_b = zeros(length(muGrid)+1,1);
fullInsurance_C_b = zeros(length(muGrid)+1,1);

parfor m=1:length(muGrid)+1
    %--------------------------------------------------------------------------
    % Parameter definition
    %--------------------------------------------------------------------------
    StaticParams = struct;
    if m==length(muGrid)+1
        mu=1;
    else
        mu = muGrid(m);
    end
    StaticParams = setModelParams(StaticParams,1,mu);
    % Use either the Proximal point algorithm (method = 0) or policy function
    % iteration (method = 1)
    StaticParams.methodEqPFI = 1; 
    
    %--------------------------------------------------------------------------
    % Initialization
    %--------------------------------------------------------------------------g
    n = length(StaticParams.kGrid_pol); 
    N = length(StaticParams.KGrid); 
    
    %--------------------------------------------------------------------------
    % Computing the solution when the aggregate state is fixed to 0 (bad state)
    %--------------------------------------------------------------------------
    c_b = repmat(log(StaticParams.kGrid_pol-StaticParams.k_min+1)+0.2,[2,1,N]);
    k_prime_b = max(StaticParams.k_min,StaticParams.wealth(0,StaticParams.KGrid,StaticParams.kGrid_pol)-c_b);
    
    [c_b,k_prime_b,diff_b]...
    = OptPol_NoAggShock_PFI(false,c_b,k_prime_b,StaticParams);
    
    % k_prime_b = k_prime_b.*(k_prime_b>StaticParams.k_min);

    %--------------------------------------------------------------------------
    % Computing the solution when the aggregate state is fixed to 1 (good state)
    %--------------------------------------------------------------------------
    if m==length(muGrid)+1
        StaticParams.mu = StaticParams.er(2)/StaticParams.er(1);
    end
    c_g = repmat(log(StaticParams.kGrid_pol-StaticParams.k_min+1)+0.2,[2,1,N]);
    k_prime_g = max(StaticParams.k_min,StaticParams.wealth(1,StaticParams.KGrid,StaticParams.kGrid_pol)-c_g);
    
    [c_g,k_prime_g,diff_g]...
    = OptPol_NoAggShock_PFI(true,c_g,k_prime_g,StaticParams);
    
    % k_prime_g = k_prime_g.*(k_prime_g>eps);
      
    %--------------------------------------------------------------------------
    % Computing the stationary distributions for the cases without aggregate
    % shocks
    %--------------------------------------------------------------------------
    [Yz,Yk] = ndgrid(1:2,StaticParams.kGrid_pol,1);
    k_prime_approx_b = @(K) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,k_prime_b,Yz,Yk,K.*ones(size(Yz)));
    k_prime_approx_g = @(K) interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,k_prime_g,Yz,Yk,K.*ones(size(Yz)));
    p_b = StaticParams.p(1:2)*2;
    p_g = StaticParams.p(3:4)*2;
    warning('off','MATLAB:eigs:SigmaNearExactEig')
    agK_stationary_b = fzero(@(K) K-calcStationaryDistr(k_prime_approx_b(K),StaticParams.P_b,p_b,StaticParams),[34,42]);
    agK_stationary_g = fzero(@(K) K-calcStationaryDistr(k_prime_approx_g(K),StaticParams.P_g,p_g,StaticParams),[34,42]);
    if m==length(muGrid)+1
        fullInsurance_K_b(m) = agK_stationary_b;
        fullInsurance_C_b(m) = interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,c_b,1,agK_stationary_b,agK_stationary_b);
        fullInsurance_K_g(m) = agK_stationary_g;
        fullInsurance_C_g(m) = interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,c_g,2,agK_stationary_g,agK_stationary_g);
    else
        [~,pdf_stationary_b] = calcStationaryDistr(k_prime_approx_b(agK_stationary_b),StaticParams.P_b,p_b,StaticParams);
        pdf_b = p_b'*pdf_stationary_b; 
        
        cdf_stationary_b = cummax(min(1,cumsum(pdf_b,2)),2);
        cdf_stationary_b(:,end) = 1;
        [cdf_k_unique,idx] = unique(cdf_stationary_b,'last');
        grid_unique = StaticParams.kGrid(idx);
        if cdf_k_unique(1)>0
            cdf_k_unique = [0,cdf_k_unique];
            grid_unique = [grid_unique(1)-1e-10,grid_unique];
        end
        quantiles_b(m,:) = interp1(cdf_k_unique,grid_unique,quantileGrid);
        
        cdf_stationary_b = cummax(min(1,cumsum(pdf_stationary_b,2)),2);
        cdf_stationary_b(:,end) = 1;
        [Yz,Yk] = ndgrid(1:2,StaticParams.kGrid,1);
        c_prime = max(0,interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,c_b,Yz,Yk,agK_stationary_b.*ones(size(Yz))));
        cdf_c_b = zeros(2,length(cdf_c_grid));
        for i=1:2
            [c_p,idx] = unique(c_prime(i,:),'last');
            cdf_c_b(i,:) = interp1([min(0,min(c_p)-1e-4),min(c_p)-1e-10,c_p,10],[0,0,cdf_stationary_b(i,idx),1],cdf_c_grid);
        end
        pdf_c_b(m,:) = p_b'*[zeros(2,1),cdf_c_b(:,2:end)-cdf_c_b(:,1:end-1)]; 
        pdf_k_b(m,:) = pdf_b;


        [~,pdf_stationary_g] = calcStationaryDistr(k_prime_approx_g(agK_stationary_g),StaticParams.P_g,p_g,StaticParams);
        pdf_g = p_g'*pdf_stationary_g; 
        
        cdf_stationary_g = cummax(min(1,cumsum(pdf_g,2)),2);
        cdf_stationary_g(:,end) = 1;
        [cdf_k_unique,idx] = unique(cdf_stationary_g,'last');
        grid_unique = StaticParams.kGrid(idx);
        if cdf_k_unique(1)>0
            cdf_k_unique = [0,cdf_k_unique];
            grid_unique = [grid_unique(1)-1e-10,grid_unique];
        end
        quantiles_g(m,:) = interp1(cdf_k_unique,grid_unique,quantileGrid);
        
        cdf_stationary_g = cummax(min(1,cumsum(pdf_stationary_g,2)),2);
        cdf_stationary_g(:,end) = 1;
        [Yz,Yk] = ndgrid(1:2,StaticParams.kGrid,1);
        c_prime = max(0,interpn(1:2,StaticParams.kGrid_pol,StaticParams.KGrid,c_g,Yz,Yk,agK_stationary_g.*ones(size(Yz))));
        cdf_c_g = zeros(2,length(cdf_c_grid));
        for i=1:2
            [c_p,idx] = unique(c_prime(i,:),'last');
            cdf_c_g(i,:) = interp1([min(0,min(c_p)-1e-4),min(c_p)-1e-10,c_p,10],[0,0,cdf_stationary_g(i,idx),1],cdf_c_grid);
        end
        pdf_c_g(m,:) = p_g'*[zeros(2,1),cdf_c_g(:,2:end)-cdf_c_g(:,1:end-1)]; 
        pdf_k_g(m,:) = pdf_g;
    end
end

compareSol.muGrid = muGrid;
compareSol.pdf_b = pdf_k_b;
compareSol.mean_b = pdf_k_b*StaticParams.kGrid';
compareSol.var_b = sum(pdf_k_b.*(repmat(StaticParams.kGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_b,[1,length(StaticParams.kGrid)])).^2,2);
compareSol.skew_b = sum(pdf_k_b.*(repmat(StaticParams.kGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_b,[1,length(StaticParams.kGrid)])).^3,2)...
                  ./max(eps,compareSol.var_b.^(3/2));
compareSol.kurt_b = sum(pdf_k_b.*(repmat(StaticParams.kGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_b,[1,length(StaticParams.kGrid)])).^4,2)...
                  ./max(eps,compareSol.var_b.^2);
compareSol.quantileGrid = quantileGrid;
compareSol.quantiles_b = quantiles_b;
compareSol.pdf_g = pdf_k_g;
compareSol.mean_g = pdf_k_g*StaticParams.kGrid';
compareSol.var_g = sum(pdf_k_g.*(repmat(StaticParams.kGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_g,[1,length(StaticParams.kGrid)])).^2,2);
compareSol.skew_g = sum(pdf_k_g.*(repmat(StaticParams.kGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_g,[1,length(StaticParams.kGrid)])).^3,2)...
                  ./max(eps,compareSol.var_g.^(3/2));
compareSol.kurt_g = sum(pdf_k_g.*(repmat(StaticParams.kGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_g,[1,length(StaticParams.kGrid)])).^4,2)...
                  ./max(eps,compareSol.var_g.^2);
compareSol.quantiles_g = quantiles_g;

compareSol.cGrid = cdf_c_grid;
compareSol.pdf_c_b = pdf_c_b;
compareSol.mean_c_b = pdf_c_b*compareSol.cGrid';
compareSol.var_c_b = sum(pdf_c_b.*(repmat(compareSol.cGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_c_b,[1,length(compareSol.cGrid)])).^2,2);
compareSol.skew_c_b = sum(pdf_c_b.*(repmat(compareSol.cGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_c_b,[1,length(compareSol.cGrid)])).^3,2)...
                  ./max(eps,compareSol.var_c_b.^(3/2));
compareSol.kurt_c_b = sum(pdf_c_b.*(repmat(compareSol.cGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_c_b,[1,length(compareSol.cGrid)])).^4,2)...
                  ./max(eps,compareSol.var_c_b.^2);
compareSol.fullInsurance_K_b = fullInsurance_K_b(end);
compareSol.fullInsurance_C_b = fullInsurance_C_b(end);
compareSol.pdf_c_g = pdf_c_g;
compareSol.mean_c_g = pdf_c_g*compareSol.cGrid';
compareSol.var_c_g = sum(pdf_c_g.*(repmat(compareSol.cGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_c_g,[1,length(compareSol.cGrid)])).^2,2);
compareSol.skew_c_g = sum(pdf_c_g.*(repmat(compareSol.cGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_c_g,[1,length(compareSol.cGrid)])).^3,2)...
                  ./max(eps,compareSol.var_c_g.^(3/2));
compareSol.kurt_c_g = sum(pdf_c_g.*(repmat(compareSol.cGrid,[length(muGrid),1])...
                     -repmat(compareSol.mean_c_g,[1,length(compareSol.cGrid)])).^4,2)...
                  ./max(eps,compareSol.var_c_g.^2);
compareSol.fullInsurance_K_g = fullInsurance_K_g(end);
compareSol.fullInsurance_C_g = fullInsurance_C_g(end);
save(strcat(folder,'/compareSol.mat'),'-struct','compareSol');