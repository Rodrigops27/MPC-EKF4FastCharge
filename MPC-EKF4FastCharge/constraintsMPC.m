function [M, gamma] = constraintsMPC(x_aug_k, mpcData)
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

    Nc = mpcData.Nc;
    M = []; gamma = [];

    useCurrent = numel(mpcData.const.constraints)>=1 && mpcData.const.constraints(1)==1;
    useVoltage = numel(mpcData.const.constraints)>=2 && mpcData.const.constraints(2)==1;
    useEta     = numel(mpcData.const.constraints)>=3 && mpcData.const.constraints(3)==1;

    %% (1) Current (absolute u) + rate (Δu) limits
    if useCurrent
        if isfield(mpcData,'pred') && isfield(mpcData.pred,'u') && isfield(mpcData.pred.u,'Cu')
            Cu = mpcData.pred.u.Cu;
        else
            Cu = tril(ones(Nc));  % accumulator
        end
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
        Phi_v = mpcData.pred.v.Phi;   G_v = mpcData.pred.v.G;
        nrows = size(Phi_v,1);        % equals q*Np (no need for q)
        v_max = mpcData.const.v_max;  v_min = mpcData.const.v_min;

        % Broadcast scalars; validate vectors
        v_max = v_max(:);  if isscalar(v_max), v_max = repmat(v_max, nrows,1); end
        v_min = v_min(:);  if isscalar(v_min), v_min = repmat(v_min, nrows,1); end
        assert(numel(v_max)==nrows && numel(v_min)==nrows, ...
               'Voltage limits must be scalar or length(size(Phi_v,1)).');

        rhs_v = Phi_v * x_aug_k;

        Mv     = [ G_v;            -G_v            ];
        gammav = [ v_max - rhs_v;
                  -(v_min - rhs_v)                  ];

        M     = [M;     Mv];
        gamma = [gamma; gammav];
    end

    %% (3) Side-reaction overpotential (η/ϕ) lower bound
    if useEta
        Phi_e = mpcData.pred.eta.Phi;   G_e = mpcData.pred.eta.G;
        nrows = size(Phi_e,1);
        eta_min = mpcData.const.phise_min;

        eta_min = eta_min(:);
        if isscalar(eta_min), eta_min = repmat(eta_min, nrows,1); end
        assert(numel(eta_min)==nrows, ...
               'phise_min must be scalar or length(size(Phi_eta,1)).');

        rhs_e = Phi_e * x_aug_k;

        Meta     = -G_e;
        gammaeta = -eta_min + rhs_e;

        M     = [M;     Meta];
        gamma = [gamma; gammaeta];
    end
end
