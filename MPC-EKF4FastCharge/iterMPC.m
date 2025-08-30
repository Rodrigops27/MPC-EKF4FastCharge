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
Sigma  = mpcData.Sigma;  % accumulator (often used in cost/constraints)
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
% y[k+i|k] = C*A^i x_k + sum_{j=0}^{i-1} C*A^{i-1-j}*B * Δu[k+1+j]  + D*u_k
[Phi_soc, G_soc,Aug_matrices] = predMat(Am,Bm,Csoc,Dsoc,Np,Nc);
mpcData.Aug_matrices_soc = Aug_matrices;

% Voltage prediction: pass Ctilde=[Cv Dv], D=0
[Phi_v, G_v,Aug_matrices] = predMat(cellState.MPC.A, cellState.MPC.B, ...
    cellState.MPC.Cv, 0, Np, Nc);
mpcData.Aug_matrices_v = Aug_matrices;

% η/ϕ prediction: pass Ctilde=[Cphi Dphi], D=0
[Phi_phi, G_phi, Aug_matrices] = predMat(cellState.MPC.A, cellState.MPC.B, ...
    cellState.MPC.Cphi, 0, Np, Nc);
mpcData.Aug_matrices_phi = Aug_matrices;

% Cache for constraints builder
mpcData.pred.soc.Phi = Phi_soc; mpcData.pred.soc.G   = G_soc;
mpcData.pred.v.Phi   = Phi_v;    mpcData.pred.v.G   = G_v;
mpcData.pred.eta.Phi = Phi_phi;  mpcData.pred.eta.G = G_phi;
mpcData.pred.u.Cu    = tril(ones(Nc));

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

[M,gamma] = constraintsMPC(dx, mpcData);

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
% mpcData.SOCk_1 = cellState.MPC.prior.SOC;
mpcData.DUk_1 = DU;

end