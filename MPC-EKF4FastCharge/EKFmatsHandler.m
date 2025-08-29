function cellState = EKFmatsHandler(ekfData, Xind, zk, Tk)
% Build Δi-augmented plant and output rows for MPC from EKF exports
% Inputs:
%   ekfData : EKF structure from iterEKF
%   Xind    : indices/weights of the 4 neighbors (from iterEKF)
%   zk      : latest EKF output vector (so we can compute Rct terms)
%   Tk      : temperature [K] used this step

% ----- weights & sizes -----
w = Xind.gamma(:);
n = ekfData.n;                 % # transient states
ind = ekfData.ind;             % indices for rows we need (Ifdl0, Ifdl3, …)
cellData = ekfData.cellData;   % ROM parameters
F = cellData.const.F;  R = cellData.const.R;

% ----- blend A (stored as diagonal vectors), then form matrix -----
A1 = ekfData.M(Xind.theT(1),Xind.theZ(1)).A;
A2 = ekfData.M(Xind.theT(2),Xind.theZ(2)).A;
A3 = ekfData.M(Xind.theT(3),Xind.theZ(3)).A;
A4 = ekfData.M(Xind.theT(4),Xind.theZ(4)).A;
Ablend_vec = w(1)*A1 + w(2)*A2 + w(3)*A3 + w(4)*A4;     % n×1
A = diag(Ablend_vec);

% Transient input gain is all ones in this ROM (see initKF), so:
B = ones(n,1);

% Δi augmentation: χ=[x; u], Δi is the input
A_aug = [A, B; zeros(1,n), 1];
B_aug = [zeros(n,1); 1];

% ----- blend C and D -----
C1 = ekfData.M(Xind.theT(1),Xind.theZ(1)).C;   D1 = ekfData.M(Xind.theT(1),Xind.theZ(1)).D;
C2 = ekfData.M(Xind.theT(2),Xind.theZ(2)).C;   D2 = ekfData.M(Xind.theT(2),Xind.theZ(2)).D;
C3 = ekfData.M(Xind.theT(3),Xind.theZ(3)).C;   D3 = ekfData.M(Xind.theT(3),Xind.theZ(3)).D;
C4 = ekfData.M(Xind.theT(4),Xind.theZ(4)).C;   D4 = ekfData.M(Xind.theT(4),Xind.theZ(4)).D;

Cb = w(1)*C1 + w(2)*C2 + w(3)*C3 + w(4)*C4;
Db = w(1)*D1 + w(2)*D2 + w(3)*D3 + w(4)*D4;

% ----- voltage row: replicate getChatV composition but include D -----
% film resistances at current SOC
x0SOC = ekfData.SOC0 - ekfData.x0*(ekfData.Ts/(3600*ekfData.Q));
SOCn  = cellData.function.neg.soc(x0SOC,Tk);
SOCp  = cellData.function.pos.soc(x0SOC,Tk);
Rfn   = cellData.function.neg.Rf(SOCn,Tk);
Rfp   = cellData.function.pos.Rf(SOCp,Tk);

% charge-transfer resistances using zk values at collectors
i0n = cellData.function.neg.k0(SOCn,Tk) ...
    * sqrt( zk(ind.Thetae0) * (1 - zk(ind.Thetass0)) * zk(ind.Thetass0) );
i0p = cellData.function.pos.k0(SOCp,Tk) ...
    * sqrt( zk(ind.Thetae3) * (1 - zk(ind.Thetass3)) * zk(ind.Thetass3) );
Rctn = R*Tk/(F*i0n);
Rctp = R*Tk/(F*i0p);

Cv =  Rfp*Cb(ind.Ifdl3,:) - Rfn*Cb(ind.Ifdl0,:) ...
    + Rctp*Cb(ind.If3,:)   - Rctn*Cb(ind.If0,:) ...
    +        Cb(ind.Phie(end),:);
Dv =  Rfp*Db(ind.Ifdl3,:) - Rfn*Db(ind.Ifdl0,:) ...
    + Rctp*Db(ind.If3,:)   - Rctn*Db(ind.If0,:) ...
    +        Db(ind.Phie(end),:);

% ----- η (side-reaction overpotential) row at negative collector -----
Cphi = Cb(ind.Phise0,:);    Dphi = Db(ind.Phise0,:);

% ----- SOC row: z = SOC0 + r * u, with r = -Ts/(3600 Q) -----
r = -ekfData.Ts/(3600*ekfData.Q);
Csoc = [zeros(1,n), r];     Dsoc = 0;

% ----- pack for MPC -----
cellState.MPC.A  = A_aug;
cellState.MPC.B  = B_aug;
cellState.MPC.Csoc = Csoc;  cellState.MPC.Dsoc = Dsoc;
cellState.MPC.Cv   = [Cv,  Dv];   % i.e., C̃_v
cellState.MPC.Dv   = 0;
cellState.MPC.Cphi = [Cphi, Dphi];% i.e., C̃_η
cellState.MPC.Dphi = 0;

% state for RHS of constraints (Φ*x_k terms):
% weighted-average transient state, and last current from mpcData later
xblend = zeros(n,1);
for k = 1:4
    xblend = xblend + w(k)*ekfData.M(Xind.theT(k),Xind.theZ(k)).xhat;
end
cellState.MPC.x_aug_k = [];           % fill as [xblend; mpcData.uk_1] before constraints
end
