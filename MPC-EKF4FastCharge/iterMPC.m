function [uk, mpcData] = iterMPC(xk,cellState,mpcData)
% ITERMPC This function computes one iteration of the MPC charging protocol
%  - Build SOC predictions (augmented Δu form)
%  - Solve unconstrained LS-QP
%  - Ask constraintsMPC for M,gamma (it will build V/eta predictions)
%  - Run Hildreth if needed
%
% Inputs:
%   xk          : present state vector
%   cellState   : structure, contains cell information from simStep
%   mpcData     : structure, contains MPC information
% Outputs:
%   uk          : optimal input current
%   mpcData     : structure, contains updated MPC information

% Load MPC data
Np     = mpcData.Np;
Nc     = mpcData.Nc;
uk_1   = mpcData.uk_1;
Ref    = mpcData.ref*ones(Np,1)/100; % targetSOC

%% Obtain matrices from Kalman Filter
% linMatrices = linMat(xk_1,mpcData);
Am     = cellState.MPC.A;
Bm     = cellState.MPC.B;
Csoc   = cellState.MPC.Csoc;
Dsoc   = cellState.MPC.Dsoc;


% Build augmented state vector
% dx = [(xk - xk_1); mpcData.SOCk_1];
dx = [xk(:); uk_1];   % uk_1 = last applied current

% Y = Phi*x + G*U
% Compute SOC prediction matrices
[Phi_soc, G_soc,Aug_matrices] = predMat(Am,Bm,Csoc,Dsoc,Np,Nc);
mpcData.Aug_matrices_soc = Aug_matrices;

% LS-QP unconstrained solution:
F = -2*G_soc'*(Ref - Phi_soc*dx);

% Adaptive control weighting
[~,S,~] = svd(G_soc'*G_soc);
[m,n] = size(S);
Ru = (norm(F,2)/(2*mpcData.const.du_max*sqrt(Nc)))-(S(m,n));
Ru = Ru*eye(mpcData.Nc,mpcData.Nc);
mpcData.Ru = Ru;

E = 2*(G_soc'*G_soc + Ru);
DU = -E\F; % Unconstrained Future inputs

[M,gamma] = constraintsMPC(dx, cellState, mpcData);

%% Check if constraints are violated; if so, run Hildreth
if sum(M*DU - gamma > 0) > 0
    [DU,lambda,nexec] = hildreth(E,F,M,gamma,mpcData.lambda,mpcData.maxHild);
    mpcData.lambda = lambda;
    mpcData.nexec = nexec;
    % disp(nexec);
else
    mpcData.nexec = 0;
end

du = DU(1);
uk = du + mpcData.uk_1;

mpcData.uk_1 = uk;
mpcData.xk_1 = xk;
% mpcData.SOCk_1 = cellState.MPC.prior.SOC;
mpcData.DUk_1 = DU;

end