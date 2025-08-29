function mpcData = initMPC(SOC0,Np,Nc,targetSOC,Nsim,ROM,constraints)
% INITMPC  Initialize MPC controller data structure.
%
% Inputs
%   SOC0       : initial SOC (percent or fraction; be consistent with the rest of the code)
%   Np         : prediction horizon (integer)
%   Nc         : control horizon (integer)
%   targetSOC  : terminal SOC reference (same units as SOC0)
%
%                .Nsim          : total sim steps (optional; stored for convenience)
%
% Outputs
%   mpcData : controller structure with fields
%             .Np, .Nc, .Ts, .ref, .Sigma
%             .Q, .Ru
%             .const.constraints, .const.(limits...)
%             .maxHild, .lambda
%             .uk_1, .SOCk_1, .xk_1, .DUk_1
%             .pred.u.Cu (and step-built .pred.soc/.pred.v/.pred.eta)

% Establish MPC tuning parameters:
mpcData.Np = Np;            % Prediction horizon
mpcData.Nc = Nc;            % Control horizon
mpcData.Ts = 1;             % Sampling interval [sec]
mpcData.ref = targetSOC;    % Target value for final SOC
mpcData.Sigma  = tril(ones(Nc)); % accumulator (often used in cost/constraints)
mpcData.Nsim = Nsim;        % Number of simulation points

% ---- Weights
mpcData.Q      = Q;
mpcData.Ru     = eye(Nc); % Example

% ---- Solver defaults / warm starts
mpcData.maxHild = 100;      % Maximum iterations for Hildreth
mpcData.lambda = [];
mpcData.uk_1 = 0;           % Initialization value
mpcData.SOCk_1 = 0;        
mpcData.DUk_1 = 0;
mpcData.SOC0 = SOC0;        % Initialization value
% mpcData.SOCk_1  = SOC0;                % last SOC

% Crate = 10; % Max Crate TBD
Iapp_max = -ROM.cellData.function.const.Q(0.5)*Crate; % Max current
Iapp_min = 0;

% Establish MPC constraints
mpcData.const.constraints = constraints; % Enter desired conditions
if constraints(1) == 1
    mpcData.const.u_max = Iapp_min;       % Minimum applied current
    mpcData.const.u_min = Iapp_max;       % Maximum applied current
    mpcData.const.du_max = 50;            % Maximum increment
    mpcData.const.du_min = -50;           % Minimum increment
end
if constraints(2) == 1
    mpcData.const.v_min = V_min;
    mpcData.const.v_max = V_max;
end

if constraints(3) == 1
    mpcData.const.phise_min = Phise_min;
end
