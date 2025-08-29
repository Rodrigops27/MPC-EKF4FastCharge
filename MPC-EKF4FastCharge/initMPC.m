function mpcData = initMPC(SOC0,Np,Nc,targetSOC,Nsim,ROM)

% Establish MPC tuning parameters:
mpcData.Np = Np;            % Prediction horizon
mpcData.Nc = Nc;            % Control horizon
mpcData.ref = targetSOC;    % Target value for final SOC

mpcData.Sigma = tril(ones(mpcData.Nc,mpcData.Nc));
mpcData.Nsim = Nsim;        % Number of simulation points

mpcData.uk_1 = 0;           % Initialization value
mpcData.SOCk_1 = 0;        
mpcData.DUk_1 = 0;

mpcData.maxHild = 100;      % Maximum iterations for Hildreth
mpcData.lambda = [];       
mpcData.SOC0 = SOC0;        % Initialization value
mpcData.Ts = 1;             % Sampling interval [sec]

Crate = 10; % Max Crate TBD
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
