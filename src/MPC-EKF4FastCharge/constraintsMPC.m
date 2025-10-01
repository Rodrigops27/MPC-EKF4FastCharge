 function [M, gamma] = constraintsMPC(x_aug_k, cellState, mpcData)
% Build stacked inequality constraints M*DU <= gamma
% Inputs:
%   x_aug_k : augmented state at time k  [ (nx+m) x 1 ]  (i.e., [x_k; u_k])
%   mpcData : struct with fields:
%       .Nc, .Np, .uk_1
%       .const.constraints = [useCurrent useVoltage useEta]
%       .const.u_min/u_max, .du_min/du_max, .v_min/.v_max, .phise_min
%       .pred.v.{Phi,G}, .pred.eta.{Phi,G}  (when enabled)
%       .pred.u.Cu (optional; defaults to tril(ones(Nc)))
%
% Output:
%   M, gamma : inequality matrices for Hildreth / QP (M*DU <= gamma)

    Nc = mpcData.Nc; Np = mpcData.Np;
    M = []; gamma = [];

    useCurrent = numel(mpcData.const.constraints)>=1 && mpcData.const.constraints(1)==1;
    useVoltage = numel(mpcData.const.constraints)>=2 && mpcData.const.constraints(2)==1;
    useEta     = numel(mpcData.const.constraints)>=3 && mpcData.const.constraints(3)==1;

    % (1) Current (absolute u) + rate (Δu) limits
    if useCurrent
        Cu = tril(ones(Nc));  % accumulator
        uk    = mpcData.uk_1;
        u_max = mpcData.const.u_max;
        u_min = mpcData.const.u_min;

        Mi     = [ Cu;            -Cu            ];
        gammai = [ (u_max - uk)*ones(Nc,1);
                  -(u_min - uk)*ones(Nc,1)      ];

        du_max = mpcData.const.du_max;
        du_min = mpcData.const.du_min;

        Mdu     = [ eye(Nc);      -eye(Nc)       ];
        gammadu = [ du_max*ones(Nc,1);
                   -du_min*ones(Nc,1)           ];

        M     = [M;     Mi;    Mdu];
        gamma = [gamma; gammai; gammadu];
    end

    %% (2) Voltage bounds
    if useVoltage
        % Voltage prediction: pass Ctilde=[Cv Dv], D=0
        [Phi_v, G_v] = predMat(cellState.MPC.A, cellState.MPC.B, ...
            cellState.MPC.Cv, cellState.MPC.Dv, Np, Nc);
        nrows = size(Phi_v,1);        % equals q*Np (no need for q)
        v_max = mpcData.const.v_max;
        
        % Broadcast scalars; validate vectors
        v_max = v_max(:);  if isscalar(v_max), v_max = repmat(v_max, nrows,1); end

        % rhs_v = Phi_v * x_aug_k + bv;
        rhs_v = Phi_v * x_aug_k + cellState.MPC.bv*ones(size(Phi_v,1),1);

        Mv     = G_v;
        gammav = v_max - rhs_v;

        M     = [M;     Mv];
        gamma = [gamma; gammav];
    end

    %% (3) Side-reaction overpotential (η/ϕ) lower bound
    if useEta
        % η/ϕ prediction: pass Ctilde=[Cphi Dphi], D=0
        [Phi_e, G_e] = predMat(cellState.MPC.A, ...
            cellState.MPC.B, cellState.MPC.Cphi, cellState.MPC.Dphi, Np, Nc);
        nrows = size(Phi_e,1);
        eta_min = mpcData.const.phise_min;

        eta_min = eta_min(:);
        if isscalar(eta_min), eta_min = repmat(eta_min, nrows,1); end
        assert(numel(eta_min)==nrows, ...
               'phise_min must be scalar or length(size(Phi_eta,1)).');

        rhs_e = Phi_e * x_aug_k + cellState.MPC.bphi*ones(size(Phi_e,1),1);

        Meta     = -G_e;
        gammaeta = -eta_min + rhs_e;

        M     = [M;     Meta];
        gamma = [gamma; gammaeta];
    end
    
    
    % --- SOC upper bound (prevents overshoot) ---
    if isfield(mpcData.const,'z_max')
        Phi_s = mpcData.Phi_soc;   % q*Np x (nx+m) with q=1 for SOC
        G_s   = mpcData.G_soc;     % q*Nc wide (Nc cols)
        rhs_s  = Phi_s * x_aug_k + mpcData.SOCk_1*ones(size(Phi_s,1),1);  % add baseline SOC_k

        zmax = mpcData.const.z_max;
        if isfield(mpcData.const,'z_tol'), zmax = zmax + mpcData.const.z_tol; end
        zmax_vec = zmax * ones(size(rhs_s));
        gam_soc_u =  zmax_vec - rhs_s;

        M     = [M;     G_s];
        gamma = [gamma; gam_soc_u];
    end

    % --- Optional: terminal lower bound (ensure reaching target) ---
    if isfield(mpcData.const,'z_target_terminal') && mpcData.const.z_target_terminal
        Phi_s = mpcData.Phi_soc;  G_s   = mpcData.G_soc;
        rN    = Phi_s(end,:)*x_aug_k;
        zt    = mpcData.const.z_max;   % reuse target as terminal min
        Mterm   = -G_s(end,:);                 % 1 x Nc
        gammaterm = -zt + rN;                  % scalar
        M     = [M;     Mterm];
        gamma = [gamma; gammaterm];
    end
end
