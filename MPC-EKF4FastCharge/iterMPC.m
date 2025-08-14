function [uk, mpcData] = iterMPC(xk,xk_1,cellState,mpcData)
% ITERMPC This function computes one iteration of the MPC charging protocol
% Inputs:
%   xk          : present state vector
%   xk_1        : previous timestep state vector
%   cellState   : structure, contains cell information from simStep
%   mpcData     : structure, contains MPC information
% Outputs:
%   uk          : optimal input current
%   mpcData     : structure, contains updated MPC information

% Load MPC data
Np     = mpcData.Np;
Nc     = mpcData.Nc;
uk_1   = mpcData.uk_1;
Sigma  = mpcData.Sigma;
Ref    = mpcData.ref*ones(Np,1)/100;

%% Obtain matrices from Kalman Filter
Am     = cellState.MPC.A;
Bm     = cellState.MPC.B;
Csoc   = cellState.MPC.Csoc;
Dsoc   = cellState.MPC.Dsoc;

% Compute SOC prediction matrices
[Phi_soc, G_soc] = predMat(Am,Bm,Csoc,Dsoc,Np,Nc);

% Build augmented state vector
dx = xk - xk_1;
SOCk_1 = cellState.MPC.SOCk_1;

%% Define values needed by Hildreth
E = -2*G_soc'*(Ref - Phi_soc*dx);
[~,S,~] = svd(G_soc'*G_soc);
[m,n] = size(S);
F = (norm(F,2)/(2*mpcData.const.du_max*sqrt(Nc)))-(S(m,n));
Ru = eye(mpcData.Nc,mpcData.Nc);
mpcData.Ru = Ru;

E = 2*(Csoc'*G_soc + Ru);
DU = -E\F;

[M,gamma] = constraintsMPC(dx,cellState,mpcData);

%% Check if constraints are violated; if so, run Hildreth
if sum(M*DU - gamma > 0) > 0
    [DU,lambda,nexec] = hildreth(E,F,M,gamma,mpcData.lambda,mpcData.maxHild);
    mpcData.lambda = lambda;
    mpcData.nexec = nexec;
else
    mpcData.nexec = 0;
end

du = DU(1);
uk = du + mpcData.uk_1;

mpcData.uk_1 = uk;
mpcData.xk_1 = xk;
mpcData.SOCk_1 = cellState.MPC.prior.SOC;
mpcData.DUk_1 = DU;

end