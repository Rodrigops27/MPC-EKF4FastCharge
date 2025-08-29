function MPC = EKFmatsHandler(ekfData, Xind, zk, Tk)
% Build Δi-augmented plant and output rows for MPC from EKF + blending.
% Outputs filled in MPC: A,B,Csoc,Dsoc,Cv,Dv,Cphi,Dphi,x_aug_k
% Inputs:
%   ekfData : EKF structure from iterEKF
%   Xind    : indices/weights of the 4 neighbors (from iterEKF)
%   zk      : latest EKF output vector (so we can compute Rct terms)
%   Tk      : temperature [K] used this step

    w = Xind.gamma(:);

    % ----- Blend A (diagonal stored as vectors), then form matrix -----
    A1 = ekfData.M(Xind.theT(1),Xind.theZ(1)).A;
    A2 = ekfData.M(Xind.theT(2),Xind.theZ(2)).A;
    A3 = ekfData.M(Xind.theT(3),Xind.theZ(3)).A;
    A4 = ekfData.M(Xind.theT(4),Xind.theZ(4)).A;

    Ablend_vec = w(1)*A1 + w(2)*A2 + w(3)*A3 + w(4)*A4;
    n = numel(Ablend_vec);
    A = diag(Ablend_vec);

    % Transient input gain (ROM design): ones
    B = ones(n,1);

    % Δi augmentation
    A_aug = [A, B; zeros(1,n), 1];
    B_aug = [zeros(n,1); 1];

    % ----- Blend C and D -----
    C1 = ekfData.M(Xind.theT(1),Xind.theZ(1)).C;  D1 = ekfData.M(Xind.theT(1),Xind.theZ(1)).D;
    C2 = ekfData.M(Xind.theT(2),Xind.theZ(2)).C;  D2 = ekfData.M(Xind.theT(2),Xind.theZ(2)).D;
    C3 = ekfData.M(Xind.theT(3),Xind.theZ(3)).C;  D3 = ekfData.M(Xind.theT(3),Xind.theZ(3)).D;
    C4 = ekfData.M(Xind.theT(4),Xind.theZ(4)).C;  D4 = ekfData.M(Xind.theT(4),Xind.theZ(4)).D;

    Cb = w(1)*C1 + w(2)*C2 + w(3)*C3 + w(4)*C4;
    Db = w(1)*D1 + w(2)*D2 + w(3)*D3 + w(4)*D4;

    % ----- Voltage row: replicate getChatV composition (include D) -----
    ind = ekfData.ind;          % structure with row indices
    cellData = ekfData.cellData;
    F = cellData.const.F;  R = cellData.const.R;

    % film resistances using operating SOC
    x0SOC = ekfData.SOC0 - ekfData.x0*(ekfData.Ts/(3600*ekfData.Q));
    SOCn  = cellData.function.neg.soc(x0SOC, Tk);
    SOCp  = cellData.function.pos.soc(x0SOC, Tk);
    Rfn   = cellData.function.neg.Rf(SOCn, Tk);
    Rfp   = cellData.function.pos.Rf(SOCp, Tk);

    % charge-transfer resistances (using zk for concentrations)
    i0n = cellData.function.neg.k0(SOCn, Tk) * sqrt( zk(ind.Thetae0) * (1 - zk(ind.Thetass0)) * zk(ind.Thetass0) );
    i0p = cellData.function.pos.k0(SOCp, Tk) * sqrt( zk(ind.Thetae3) * (1 - zk(ind.Thetass3)) * zk(ind.Thetass3) );
    Rctn = R*Tk/(F*i0n);
    Rctp = R*Tk/(F*i0p);

    Cv =  Rfp*Cb(ind.Ifdl3,:) - Rfn*Cb(ind.Ifdl0,:) ...
        + Rctp*Cb(ind.If3,:)   - Rctn*Cb(ind.If0,:) ...
        +        Cb(ind.Phie(end),:);
    Dv =  Rfp*Db(ind.Ifdl3,:) - Rfn*Db(ind.Ifdl0,:) ...
        + Rctp*Db(ind.If3,:)   - Rctn*Db(ind.If0,:) ...
        +        Db(ind.Phie(end),:);

    % ----- η (side-reaction overpotential) at negative collector -----
    Cphi = Cb(ind.Phise0,:); 
    Dphi = Db(ind.Phise0,:);

    % ----- SOC row: z = SOC0 + r * u   with r = -Ts/(3600*Q)
    r = -ekfData.Ts/(3600*ekfData.Q);
    Csoc = [zeros(1,n), r];
    Dsoc = 0;

    % ----- Pack results -----
    MPC.A  = A_aug;
    MPC.B  = B_aug;

    MPC.Csoc = Csoc;                MPC.Dsoc = Dsoc;
    MPC.Cv   = [Cv,  Dv];           MPC.Dv   = 0;
    MPC.Cphi = [Cphi, Dphi];        MPC.Dphi = 0;

    % % initial augmented state for RHS terms (x_k; u_k)
    % if isfield(ekfData,'xhat')
    %     xk = ekfData.xhat(:);
    % else
    %     % fallback: weighted average of neighbor xhats if present
    %     xk = zeros(n,1);
    %     for kk=1:4
    %         if isfield(ekfData.M(Xind.theT(kk),Xind.theZ(kk)), 'xhat')
    %             xk = xk + w(kk) * ekfData.M(Xind.theT(kk),Xind.theZ(kk)).xhat(:);
    %         end
    %     end
    % end
    % uk = mpcData.uk_1;
    % MPC.x_aug_k = [xk; uk];
end
