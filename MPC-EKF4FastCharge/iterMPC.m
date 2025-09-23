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
Ref    = mpcData.ref*ones(Np,1); % targetSOC

%% Obtain matrices from Kalman Filter
Am     = cellState.MPC.A;
Bm     = cellState.MPC.B;
Csoc   = cellState.MPC.Csoc;
Dsoc   = cellState.MPC.Dsoc;

% Build augmented state vector
dx = [xk(:); uk_1];   % uk_1 = last applied current

% Y = Phi*x + G*U
% Compute SOC prediction matrices
[Phi_soc, G_soc, Aug_matrices] = predMat(Am,Bm,Csoc,Dsoc,Np,Nc);
mpcData.Phi_soc = Phi_soc; mpcData.G_soc = G_soc;
% Aug_matrices useful for quick stability analysis

% LS-QP unconstrained solution:
F = -2*G_soc'*(Ref - Phi_soc*dx);

% Adaptive control weighting
[~,S,~] = svd(G_soc'*G_soc); % lambda min estimation
[m,n] = size(S); 
Ru = (norm(F,2)/(2*mpcData.const.du_max*sqrt(Nc)))-(S(m,n));
Ru = Ru*eye(Nc,Nc);
mpcData.Ru = Ru;

E = 2*(G_soc'*G_soc + Ru);
DU = -E\F; % Unconstrained Future inputs

J_unc = (Ref - Phi_soc*dx - G_soc*DU)'*(Ref - Phi_soc*dx - G_soc*DU) + ...
    DU'*Ru*DU;

% Stability Analysis
Kmpc = E \ (G_soc'*Phi_soc);
Kmpc = Kmpc(1,:);
CL = (Aug_matrices.A - Aug_matrices.B*Kmpc);
[~,poles] = eig(CL);
poles = diag(poles);
mpcData.poles = poles;
mpcData.sv = svd(CL);
% disp(poles);

[M,gamma] = constraintsMPC(dx, cellState, mpcData);

%% Check if constraints are violated; if so, run Hildreth
if sum(M*DU - gamma > 0) > 0
    [DU,lambda,nexec] = hildreth(E,F,M,gamma,mpcData.lambda,mpcData.maxHild);
    mpcData.lambda = lambda; % TBD: faster without it?
    mpcData.nexec = nexec;

else
    mpcData.nexec = 0;
end

du = DU(1);
uk = du + mpcData.uk_1;

mpcData.uk_1 = uk;
mpcData.xk_1 = xk;
mpcData.DUk_1 = DU;

vres   = M*DU - gamma;
nviol  = sum(vres > 1e-9);

J_fin = (Ref - Phi_soc*dx - G_soc*DU)'*(Ref - Phi_soc*dx - G_soc*DU) + ...
    DU'*Ru*DU;

% Cost logging
kidx = mpcData.k;
mpcData.cost.t(kidx,1)        = (kidx-1) * mpcData.Ts;
mpcData.cost.J_uncon(kidx,1)  = J_unc;
mpcData.cost.J_final(kidx,1)  = J_fin;
mpcData.cost.norm_DU(kidx,1)  = norm(DU,2);
mpcData.cost.viol(kidx,1)     = nviol;
mpcData.cost.nexec(kidx,1)    = mpcData.nexec;
end