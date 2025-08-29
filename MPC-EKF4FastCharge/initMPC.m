function mpcData = initMPC(SOC0, Np, Nc, targetSOC, opts)
% INITMPC  Initialize MPC controller data structure.
%
% Inputs
%   SOC0       : initial SOC (percent or fraction — be consistent project-wide)
%   Np         : prediction horizon (integer)
%   Nc         : control horizon (integer)   (will be clamped to Np)
%   targetSOC  : terminal SOC reference (same units as SOC0)
%   opts       : (struct, optional) configuration
%                .Ts            : sample time [s]                          (default: 1)
%                .constraints   : [useCurrent useVoltage useEta] logical   (default: [1 1 1])
%                .limits        : struct with any of:
%                                  - u_min, u_max        (A)
%                                  - du_min, du_max      (A/sample)
%                                  - v_min, v_max        (V or Np×1)
%                                  - phise_min           (V or Np×1)
%                .Q, .Ru        : weights (scalar or matrices)             (defaults: Q=1, Ru=eye(Nc))
%                .Nsim          : total sim steps (optional; stored only for convenience)
%                .ROM, .Crate   : if provided and u_min/u_max missing, they’re used to derive current limits
%
% Output
%   mpcData : controller structure with fields
%             .Np, .Nc, .Ts, .ref, .Sigma
%             .const.constraints, .const.(limits...)
%             .maxHild, .lambda
%             .uk_1, .SOCk_1, .xk_1, .DUk_1
%             .pred.u.Cu   (accumulator; other .pred.* filled each step)

if nargin < 5, opts = struct(); end

% Basic sizes / clamps
Np = max(1, round(Np));
Nc = max(1, min(round(Nc), Np));

% ---- Basics
mpcData.Np     = Np;
mpcData.Nc     = Nc;
mpcData.Ts     = getfield_with_default(opts, 'Ts', 1);
mpcData.ref    = targetSOC;            % terminal SOC target
mpcData.Sigma  = tril(ones(Nc));       % accumulator (often used in cost/constraints)
mpcData.Nsim   = getfield_with_default(opts, 'Nsim', NaN);

% Output Weight?
% mpcData.Q      = getfield_with_default(opts, 'Q',  1);

% Solver / warm starts
mpcData.maxHild = 100;
mpcData.lambda  = [];                  % dual warm start
mpcData.DUk_1   = zeros(Nc,1);         % last optimal Δu sequence (optional)
mpcData.xk_1    = [];                  % last state (optional)
mpcData.uk_1    = 0;                   % last applied current
% mpcData.SOCk_1  = SOC0;                % last SOC
mpcData.SOC0    = SOC0;

% Establish MPC constraints
cons = getfield_with_default(opts, 'constraints', [1 1 1]);  % [current voltage eta]
cons = logical(cons(:)).';
% if numel(cons) < 3, cons(end+1:3) = false; end
mpcData.const.constraints = cons;

% Limits (merge provided limits with safe defaults)
L = struct;
if isfield(opts, 'limits'), L = opts.limits; end

% Derive Iapp_max from ROM + Crate
Imax = -opts.ROM.cellData.function.const.Q(0.5)*opts.Crate;  % magnitude [A]
L.u_min = Imax;

% Store limits
% if ~isfield(mpcData,'const'), mpcData.const = struct; end
mpcData.const = copyfields(mpcData.const, L);

% ---- Prediction-cache seeds
mpcData.pred = struct();
mpcData.pred.u.Cu = tril(ones(Nc));    % for absolute-u constraints; others built per step
end

% === helpers ===
function v = getfield_with_default(S, name, default)
if isfield(S,name), v = S.(name); else, v = default; end
end

function A = copyfields(A, B)
f = fieldnames(B);
for k = 1:numel(f), A.(f{k}) = B.(f{k}); end
end
