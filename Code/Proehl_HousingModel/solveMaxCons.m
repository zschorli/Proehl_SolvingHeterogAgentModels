% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This file finds maximum feasible consumption c* given wealth and housing.
%__________________________________________________________________________

function [A, J] = solveMaxCons(sol, wealth, uq, h_bound, SymFuns)
%   Solves: F_Cond(c) = wealth - uq * min(h_bound, h_uc(c)) - c = 0
%   where:  h_uc(c) = eta * uq^(-eps) * c
%
%   sol     : scalar current guess for c
%   wealth  : scalar total wealth
%   uq      : scalar user cost = (1+r)*q - q'
%   h_bound : scalar collateral bound
%   StaticParams : model parameters (unused here — kept for interface compat.)
%   SymFuns : precomputed handles from buildSymFuns
%
%   Uses smooth_min consistent with solveEquCond and buildSymFuns.

    % Forward pass: F_Cond(c) = wealth - uq*h(c) - c
    % SymFuns.MaxCons.F_Cond returns h = h_uc(c)
    % Argument order: (c, uq, hb)

    % mu = 1e-4;   % must match buildSymFuns
    % smooth_max = @(x) max(x, 0) + mu .* log(1 + exp(-abs(x) ./ mu));
    % smooth_min = @(x, y) x-smooth_max(x-y);

    h_val = SymFuns.MaxCons.F_Cond(sol, uq, h_bound);
    A = wealth - uq .* h_val - sol;

    if nargout > 1
        % Derivative: dF/dc = -uq * dh/dc - 1
        % SymFuns.MaxCons.dF_dc returns d/dc[-uq*h(c) - c] = -uq*dh/dc - 1
        % (wealth is constant w.r.t. c, so it drops out)
        J = SymFuns.MaxCons.dF_dc(sol, uq, h_bound);
    end
end