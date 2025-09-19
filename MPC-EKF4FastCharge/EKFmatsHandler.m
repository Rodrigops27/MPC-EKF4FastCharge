function [MPC, xhat] = EKFmatsHandler(ekfData, Xind, zk, Tk)
% EKFmatsHandler  (nearest model, ORIGINAL SS ONLY)
% Pick the closest local ROM (no blending) and return the linearized
% state-space needed by MPC/constraints *without* Δu augmentation.
%
% Inputs
%   ekfData : EKF struct with a grid of local models ekfData.M(iT,iZ)
%             Each M(iT,iZ) has fields A,B (or A as diag vector), C, D.
%   Xind    : struct with fields:
%                .theT(1x4), .theZ(1x4)  indices of the 4 neighbors
%                .gamma(1x4)             bilinear weights (sum=1)
%            We choose the corner with max gamma (nearest).
%   zk      : latest EKF output vector (needed for voltage composition)
%   Tk      : temperature at this step (same units used by ROM functions)
%
% Output
%   MPC : struct with ORIGINAL (non-augmented) SS rows:
%           .A, .B                  % nx×nx, nx×m
%           .Csoc, .Dsoc            % 1×nx, 1×m      (SOC row)
%           .Cv,   .Dv              % 1×nx, 1×m      (voltage row)
%           .Cphi, .Dphi            % 1×nx, 1×m      (side-reaction η row)
%         plus trace fields:
%           .iT, .iZ, .pickWeight

    % --------- 1) choose nearest corner by largest gamma ----------
    [~, imax] = max(Xind.gamma(:));
    iT = Xind.theT(imax);
    iZ = Xind.theZ(imax);

    mdl = ekfData.M(iT, iZ);

    % --------- 2) Setting the state space model
    xhat = [mdl.xhat; ekfData.xhat(end)];  % + Integrator state

    mdl.A = [mdl.A; 1]; % + Integrator state
    nx = numel(mdl.A); 
    A  = diag(mdl.A(:));
    B = ones(nx,1);
    C = mdl.C;
    D = mdl.D; 

    % ----- SOC row: z = SOC0 + r * u   with r = -Ts/(3600*Q)
    r = -ekfData.Ts/(3600*ekfData.Q);
    Csoc = [zeros(1,nx-1), r];
    Dsoc = 0;

    % 3b) Voltage row: v_k = Cv*x_k + Dv*u_k + b_v,k
    %     Compose from linear outputs using ROM indices if available.
    ind      = ekfData.ind;           % structure of output row indices
    cellData = ekfData.cellData;
    F        = cellData.const.F;
    R        = cellData.const.R;
    TK = Tk+273.15;

    % Film resistances at current operating SOC
    SOCavg  = ekfData.SOC0;          % use EKF's running SOC0 (project conv.)
    SOCn    = cellData.function.neg.soc(SOCavg, TK);
    SOCp    = cellData.function.pos.soc(SOCavg, TK);
    Rfn     = cellData.function.neg.Rf(SOCn, TK);
    Rfp     = cellData.function.pos.Rf(SOCp, TK);

    % Charge-transfer "resistances" via exchange-current i0 at collectors
    i0n = cellData.function.neg.k0(SOCn, TK) * ...
        sqrt( zk(ind.Thetae0) * (1 - zk(ind.Thetass0)) * zk(ind.Thetass0) );
    i0p = cellData.function.pos.k0(SOCp, TK) * ...
        sqrt( zk(ind.Thetae3) * (1 - zk(ind.Thetass3)) * zk(ind.Thetass3) );
    Rctn = R*TK/(F*i0n);
    Rctp = R*TK/(F*i0p);

    % Compose voltage row from linear outputs (Ifdl/If/Phie)
    Cv =  Rfp*C(ind.Ifdl3,:) - Rfn*C(ind.Ifdl0,:) ...
        + Rctp*C(ind.If3,:)  - Rctn*C(ind.If0,:)  ...
        +       C(ind.Phie(end),:);

    Dv =  Rfp*D(ind.Ifdl3,:) - Rfn*D(ind.Ifdl0,:) ...
        + Rctp*D(ind.If3,:)  - Rctn*D(ind.If0,:)  ...
        +       D(ind.Phie(end),:);

    % ---- Voltage bias term b_v,k ----
    Upos = cellData.function.pos.Uocp(SOCp, TK);   % U_ocp^pos(c_s,e,k)
    Uneg = cellData.function.neg.Uocp(SOCn, TK);   % U_ocp^neg(c_s,e,k(0))
    bphi = 0;  % optional electrolyte bias; set to 0 unless you have it explicitly
    bv   = (Upos - Uneg) + bphi;


    % 3c) Side-reaction overpotential row at negative collector:
    Cphi = zeros(1,nx); Dphi = 0;
    try
        Cphi = C(ind.Phise0,:);
        Dphi = D(ind.Phise0,:);
    catch
        % leave zeros if index not present
    end


    % --------- 4) Pack result (ORIGINAL SS only) ----------
    MPC = struct();
    MPC.A    = A;
    MPC.B    = B;

    MPC.Csoc = Csoc;  MPC.Dsoc = Dsoc;
    MPC.Cv   = [Cv 0];    MPC.Dv   = Dv;
    MPC.Cphi = [Cphi 0];  MPC.Dphi = Dphi;
    MPC.bv = bv;    % Voltage bias term

    % Traceability
    MPC.iT = iT;
    MPC.iZ = iZ;
    MPC.pickWeight = Xind.gamma(imax);
end
