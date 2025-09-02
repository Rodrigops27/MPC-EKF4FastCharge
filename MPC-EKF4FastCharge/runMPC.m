% runMPC: Runs an MPC to compute the optimal current for a fast charge

clear; clc;
ROMfile = 'ROM_NMC30_HRA.mat';    % generated ROM
load(ROMfile,'ROM');

% Simulation conditions
TC = 25;                          % [°C] simulation temperature
Ts = 1;                           % [s] sample time
finalTime = 1200;                 % [s]
time = (0:finalTime).';           % N samples as a column
N = numel(time);
SOC0 = 10;                        % [%] initial SOC (project-wide convention)

%% EKF data
blend   = 'OutB';
SigmaX0 = diag([ones(1,5) 2e6]);  % initial state uncertainty
SigmaW  = 1e2;                    % process/current noise
SigmaV  = 1e-3;                   % voltage noise
ekfData = initKF(1*SOC0, TC, SigmaX0, SigmaV, SigmaW, blend, ROM);
ekfData.maxWarn  = 1000;
ekfData.trueSOC0 = SOC0/100;
% Reserve storage for EKF results... including voltage and SOC (the "+2")
zkEst   = NaN(ekfData.nz+2, N);
zkBound = NaN(size(zkEst));

%% MPC data
Np = 5;            % Prediction horizon
Nc = 2;            % Control horizon
targetSOC = 95;    % [%] terminal target

% Constraint switches (current, voltage, eta/phi)
constraints = [1 0 0];

% Limits
Crate = 1;                          % fast-charge C-rate
limits.u_max     = 0;               % forbid discharge during fast charge
limits.du_min    = -50;             % A/sample
limits.du_max    =  50;
limits.v_min     = 3.0;             % V
limits.v_max     = 4.3;             % V
limits.phise_min = 0.08;            % V (80 mV)
% mpcData.const.z_max = SOC_ref/100;  % Max SOC

opts = struct('Ts', Ts, ...
              'constraints', constraints, ...
              'limits', limits, ...
              'ROM', ROM, ...
              'Crate', Crate, ...
              'Nsim', finalTime);

mpcData = initMPC(SOC0, Np, Nc, targetSOC, opts);

% Storing variables
u_applied      = NaN(N,1);          % actually applied to ROM at k
u_store        = NaN(N,1);          % MPC command computed at k (for k+1)
voltage_store  = NaN(N,1);
x_store        = NaN(1,N);          % SOC trace for plotting

% Plant initialization
uk=0;
initCfg = struct('SOC0', SOC0, 'warnOff', true);
[voltage, obs, cellState] = OB_step(uk, TC, [], ROM, initCfg);
% Seed k=1
voltage_store(1) = voltage;
u_store(1) = uk;
x_store(1,1)     = SOC0;            % or zk(end) if you run EKF at k=1


% -------------------------------
% Main time loop (k = 1..N)
% -------------------------------
for k = 1:N
    % 1) Apply current to ROM at this step
    [voltage, obs, cellState] = OB_step(uk, TC, cellState, ROM);

    % 2) EKF update
    [zk, zbk, ekfData, Xind] = iterEKF(voltage, uk, TC, ekfData);

    % 3) EKF -> MPC linearization (Δi form)
    [MPC, xhat] = EKFmatsHandler(ekfData, Xind, zk, TC);
    cellState.MPC = MPC;
    % 4) MPC inputs
    xk = xhat;              % state estimate at k
    mpcData.SOCk_1 = zk(end);

    % 5) MPC step -> command uk for NEXT step
    [uk, mpcData] = iterMPC(xk, cellState, mpcData);

    % 6) Storing
    u_store(k)       = uk;         % MPC command
    voltage_store(k) = voltage;
    x_store(1,k)     = zk(end);
    zkEst(:,k)       = zk;
    zkBound(:,k)     = zbk;

end

% Final plots (plotFigures expects: time, u_store, voltage_store, x_store, mpcData)
plotFigures;
