% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This file solves Market clearing residuals for housing (H=1) and 
% assets (A=q').
%__________________________________________________________________________

function [A, J] = solveEquCond(sol, q, r, a, aux, h_bound, basic_grid, qmin, qmax, StaticParams, Poly, SymFuns)
% sol  = [uq_b; uq_g; qpr_b; qpr_g]   (4x1)
%          uq  = (1+r)*q_cur - q'      user cost of housing
%          qpr = next-period housing price
%
% q        : (4 x nA)  current housing price, rows = [bu, be, gu, ge]
% a        : (4 x nA)  individual asset grid
% aux      : (4 x nA)  Euler equation RHS (from OptPol loop)
% h_bound  : (4 x nA)  collateral bound = a / ((1-theta)*q)
% basic_grid, StaticParams, Poly : passed through to calcAggregates
% SymFuns  : precomputed derivative handles from buildSymFuns
%
% Output:
%   A = [H_b-1; H_g-1; A_b-qpr_b; A_g-qpr_g]   (4x1)
%   J = dA/d(sol)                                 (4x4)  [only if nargout>1]

    mu = 1e-2;   % must match buildSymFuns
    smooth_max = @(x) max(x, 0) + mu .* log(1 + exp(-abs(x) ./ mu));
    smooth_min = @(x, y) x-smooth_max(x-y);

    rpr_b  = sol(1);
    rpr_g  = sol(2);
    qpr_b = sol(3);
    qpr_g = sol(4);

   % =====================================================================
    % FORWARD PASS
    % =====================================================================

    q_b = q(1, :);
    q_g = q(3, :);
    r_b = r(1, :);
    r_g = r(3, :);

    % Raw user-cost arguments (before smoothing)
    v_b = (1+rpr_b).*q_b - qpr_b;   % (1 x nA)
    v_g = (1+rpr_g).*q_g - qpr_g;

    % Smoothed user cost
    uq_b = smooth_max(v_b);
    uq_g = smooth_max(v_g);

    % OptC_mult
    OptC_b = (1 + StaticParams.eta .* uq_b.^(1 - StaticParams.eps)) ...
             .^((1 - StaticParams.eps*StaticParams.gamma) ./ (StaticParams.eps - 1));
    OptC_g = (1 + StaticParams.eta .* uq_g.^(1 - StaticParams.eps)) ...
             .^((1 - StaticParams.eps*StaticParams.gamma) ./ (StaticParams.eps - 1));

    % Consumption
    c_b = aux(1:2, :) .*(OptC_b).^(1./StaticParams.gamma);
    c_g = aux(3:4, :) .*(OptC_g).^(1./StaticParams.gamma);

    % Housing demand
    h_uc_b = StaticParams.eta .* uq_b.^(-StaticParams.eps) .* c_b;
    h_uc_g = StaticParams.eta .* uq_g.^(-StaticParams.eps) .* c_g;
    h_b = smooth_min(h_bound(1:2, :), h_uc_b);
    h_g = smooth_min(h_bound(3:4, :), h_uc_g);

    % Wealth
    wealth_b = StaticParams.wage(1:2, :) + (1+r_b) .* a(1:2, :);
    wealth_g = StaticParams.wage(3:4, :) + (1+r_g) .* a(3:4, :);

    % Asset savings
    a_b = smooth_max(wealth_b - c_b - uq_b .* h_b);
    a_g = smooth_max(wealth_g - c_g - uq_g .* h_g);

    % Interpolate and aggregate
    h_all = [h_b; h_g];
    a_all = [a_b; a_g];
    nFine = length(StaticParams.aGrid);
    aGrid_safe = min(max(StaticParams.aGrid, StaticParams.a_min), StaticParams.aGrid_pol(end));
    h_large = interp1(StaticParams.aGrid_pol, h_all', aGrid_safe)';
    a_large = interp1(StaticParams.aGrid_pol, a_all', aGrid_safe)';

    H = calcAggregates(h_large, basic_grid, StaticParams, Poly);
    A_raw = calcAggregates(a_large, basic_grid, StaticParams, Poly);
    % Smoothly cap aggregate asset/price-implied quantity A componentwise.
    A = qmin + smooth_max(A_raw - qmin) - smooth_max(A_raw - qmax);

    A = [H - 1; A - sol(3:4)];

    % =====================================================================
    % JACOBIAN
    % sol = [rpr_b; rpr_g; qpr_b; qpr_g]
    % =====================================================================
   if nargout > 1

        sigmoid = @(x) 1 ./ (1 + exp(-x ./ mu));

        % -----------------------------------------------------------------
        % Sigmoid factors
        %
        % sig_v_b = d(uq_b)/d(v_b) = sigmoid(v_b/mu)   (1 x nA)
        % sig_a_b = d(smax)/d(arg) at savings argument  (2 x nA)
        %
        % Chain rules (both rpr and qpr act only through v):
        %   d(v_b)/d(rpr_b) =  q_b    → scale sig_v by  q_b
        %   d(v_b)/d(qpr_b) = -1      → scale sig_v by -1
        % -----------------------------------------------------------------
        sig_v_b = sigmoid(v_b);                          % (1 x nA)
        sig_v_g = sigmoid(v_g);
        sig_a_b = sigmoid(wealth_b - c_b - uq_b .* h_b);  % (2 x nA)
        sig_a_g = sigmoid(wealth_g - c_g - uq_g .* h_g);

        % -----------------------------------------------------------------
        % Symbolic derivatives w.r.t. v
        % SymFuns argument order: (v, a, hb)
        % -----------------------------------------------------------------
        dh_b_dv     = SymFuns.dh_dv(    v_b, aux(1:2,:), h_bound(1:2,:));  % (2 x nA)
        da_raw_b_dv = SymFuns.da_raw_dv(v_b, aux(1:2,:), h_bound(1:2,:));
        dh_g_dv     = SymFuns.dh_dv(    v_g, aux(3:4,:), h_bound(3:4,:));
        da_raw_g_dv = SymFuns.da_raw_dv(v_g, aux(3:4,:), h_bound(3:4,:));

        % da/dv includes smax chain rule at savings
        da_b_dv = sig_a_b .* da_raw_b_dv;
        da_g_dv = sig_a_g .* da_raw_g_dv;

        % -----------------------------------------------------------------
        % Apply outer chain rule: d(v)/d(sol)
        %
        %   d/d(rpr_b):  d(v_b)/d(rpr_b) = q_b   (1 x nA, broadcasts over rows)
        %   d/d(qpr_b):  d(v_b)/d(qpr_b) = -1
        %   cross terms: v_b independent of rpr_g, qpr_g → zero
        % -----------------------------------------------------------------
        dh_b_drprb  =  sig_v_b .* q_b .* dh_b_dv;      % (2 x nA)
        dh_b_dqprb  = -sig_v_b        .* dh_b_dv;
        da_b_drprb  =  sig_v_b .* q_b .* da_b_dv;
        da_b_dqprb  = -sig_v_b        .* da_b_dv;

        dh_g_drprg  =  sig_v_g .* q_g .* dh_g_dv;
        dh_g_dqprg  = -sig_v_g        .* dh_g_dv;
        da_g_drprg  =  sig_v_g .* q_g .* da_g_dv;
        da_g_dqprg  = -sig_v_g        .* da_g_dv;

        % -----------------------------------------------------------------
        % Interpolate onto fine aGrid
        % -----------------------------------------------------------------
        zeros24   = zeros(2, nFine);
        interp_up = @(x) interp1(StaticParams.aGrid_pol, x', aGrid_safe)';

        dh_b_drprb_l = interp_up(dh_b_drprb);
        dh_g_drprg_l = interp_up(dh_g_drprg);
        dh_b_dqprb_l = interp_up(dh_b_dqprb);
        dh_g_dqprg_l = interp_up(dh_g_dqprg);
        da_b_drprb_l = interp_up(da_b_drprb);
        da_g_drprg_l = interp_up(da_g_drprg);
        da_b_dqprb_l = interp_up(da_b_dqprb);
        da_g_dqprg_l = interp_up(da_g_dqprg);

        % -----------------------------------------------------------------
        % Aggregate via PCE
        % -----------------------------------------------------------------
        dH_drprb = calcAggregates([dh_b_drprb_l; zeros24], basic_grid, StaticParams, Poly);
        dH_drprg = calcAggregates([zeros24; dh_g_drprg_l], basic_grid, StaticParams, Poly);
        dH_dqprb = calcAggregates([dh_b_dqprb_l; zeros24], basic_grid, StaticParams, Poly);
        dH_dqprg = calcAggregates([zeros24; dh_g_dqprg_l], basic_grid, StaticParams, Poly);

        dAraw_drprb   = calcAggregates([da_b_drprb_l;   zeros24], basic_grid, StaticParams, Poly);
        dAraw_drprg   = calcAggregates([zeros24; da_g_drprg_l  ], basic_grid, StaticParams, Poly);
        dAraw_dqprb = calcAggregates([da_b_dqprb_l; zeros24], basic_grid, StaticParams, Poly);
        dAraw_dqprg = calcAggregates([zeros24; da_g_dqprg_l], basic_grid, StaticParams, Poly);

        % -------------------------------------------------------------
        % Chain rule through the smooth aggregate clamp:
        %
        % A = qmin + smax(A_raw-qmin) - smax(A_raw-qmax)
        %
        % dA/dA_raw = sigmoid(A_raw-qmin) - sigmoid(A_raw-qmax).
        %
        % A_raw, dAraw_* and dAcap_dAraw are all (2 x 1), corresponding
        % to the bad and good aggregate states.
        % -------------------------------------------------------------
        dAcap_dAraw = sigmoid(A_raw - qmin) - sigmoid(A_raw - qmax);

        dA_drprb   = dAcap_dAraw .* dAraw_drprb;
        dA_drprg   = dAcap_dAraw .* dAraw_drprg;
        dA_dqprb = dAcap_dAraw .* dAraw_dqprb;
        dA_dqprg = dAcap_dAraw .* dAraw_dqprg;

        % -----------------------------------------------------------------
        % Assemble 4x4 Jacobian
        %
        % A   = [H_b - 1    ]     sol = [rpr_b]
        %       [H_g - 1    ]           [rpr_g]
        %       [A_b - qpr_b]           [qpr_b]
        %       [A_g - qpr_g]           [qpr_g]
        %
        %            d/drpr_b      d/drpr_g     d/dqpr_b        d/dqpr_g
        % dH_b :  dH_drprb(1)   dH_drprg(1)  dH_dqprb(1)    dH_dqprg(1)
        % dH_g :  dH_drprb(2)   dH_drprg(2)  dH_dqprb(2)    dH_dqprg(2)
        % dA_b :  dA_drprb(1)   dA_drprg(1)  dA_dqprb(1)-1  dA_dqprg(1)
        % dA_g :  dA_drprb(2)   dA_drprg(2)  dA_dqprb(2)    dA_dqprg(2)-1
        % -----------------------------------------------------------------
        J = [ dH_drprb(1),  dH_drprg(1),  dH_dqprb(1),    dH_dqprg(1);
              dH_drprb(2),  dH_drprg(2),  dH_dqprb(2),    dH_dqprg(2);
              dA_drprb(1),  dA_drprg(1),  dA_dqprb(1)-1,  dA_dqprg(1);
              dA_drprb(2),  dA_drprg(2),  dA_dqprb(2),    dA_dqprg(2)-1 ];
    end
end