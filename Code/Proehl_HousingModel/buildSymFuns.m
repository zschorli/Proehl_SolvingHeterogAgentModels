% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% Precompute symbolic derivatives for solveEquCond and solveMaxCons.
%__________________________________________________________________________

function SymFuns = buildSymFuns(StaticParams, mu)
%   sol = [r_b; r_g; qpr_b; qpr_g]
%
%   v        = (1+r)*q - qpr          raw user-cost argument
%   uq       = smooth_max(v)          smoothed user cost
%   c        = aux * OptC(uq)^(1/gam) consumption
%   wealth   = wage + (1+r)*a         wealth depends on r directly
%   a_raw    = wealth - c - uq*h      savings (wealth term added)
%
%   Both r and qpr enter through v:
%     d(v)/d(r)   =  q   → d(uq)/d(r)   = sigmoid(v)*q
%     d(v)/d(qpr) = -1   → d(uq)/d(qpr) = -sigmoid(v)
%
%   r also enters wealth directly:
%     d(wealth)/d(r) = a   → extra term in da_raw/dr
%
%   So we need:
%     d/dv  : through uq (same as before)
%     da_raw_da_direct : the extra +a term from wealth in da_raw/dr

eps_val = StaticParams.eps;
gam     = StaticParams.gamma;
eta     = StaticParams.eta;

sym_smax = @(x)   x + mu .* log(1 + exp(-x ./ mu));
sym_smin = @(x,y) x - sym_smax(x - y);

% =========================================================================
% PART 1 — Derivatives for solveEquCond
%
% Symbolic variables:
%   v  : (1+r)*q - qpr    raw user-cost argument
%   a  : aux              (parameter)
%   hb : h_bound          (parameter)
%
% Note: wealth = wage + (1+r)*a is constant w.r.t. v and qpr,
%       but has a direct derivative w.r.t. r via a (handled numerically).
%
% a_raw = -c(v) - uq(v)*h(v)   [wealth dropped — added back numerically]
% =========================================================================

syms v a hb real

uq    = sym_smax(v);
OptC  = (1 + eta .* uq.^(1 - eps_val)).^((1 - eps_val*gam)./(eps_val - 1));
F_c   = a .* OptC.^(1./gam);

F_h_uc = eta .* uq.^(-eps_val) .* F_c;
F_h    = sym_smin(hb, F_h_uc);

% a_raw excludes wealth (wealth's direct r-dependence handled numerically)
F_a_raw = -F_c - uq .* F_h;

dh_dv     = diff(F_h,     v);
da_raw_dv = diff(F_a_raw, v);

vars = [v, a, hb];

raw.F_h       = matlabFunction(F_h,       'Vars', vars);
raw.F_c       = matlabFunction(F_c,       'Vars', vars);
raw.dh_dv     = matlabFunction(dh_dv,     'Vars', vars);
raw.da_raw_dv = matlabFunction(da_raw_dv, 'Vars', vars);

fns = fieldnames(raw);
for i = 1:numel(fns)
    f = raw.(fns{i});
    SymFuns.(fns{i}) = @(varargin) safeEval(f, varargin{:});
end

% =========================================================================
% PART 2 — Derivatives for solveMaxCons (unchanged)
% =========================================================================
syms c_sym uq_mc hb_mc real

F_h_uc_mc = eta .* uq_mc.^(-eps_val) .* c_sym;
F_h_mc    = sym_smin(hb_mc, F_h_uc_mc);
F_Cond_mc = -uq_mc .* F_h_mc - c_sym;
dF_dc_mc  = diff(F_Cond_mc, c_sym);

vars_mc = [c_sym, uq_mc, hb_mc];

SymFuns.MaxCons.F_Cond = @(varargin) safeEval(matlabFunction(F_h_mc,   'Vars', vars_mc), varargin{:});
SymFuns.MaxCons.dF_dc  = @(varargin) safeEval(matlabFunction(dF_dc_mc, 'Vars', vars_mc), varargin{:});

% fprintf('buildSymFuns: precomputed (solveEquCond + solveMaxCons).\n');
end

function out = safeEval(f, varargin)
    out = f(varargin{:});
    out(~isfinite(out)) = 0;
end