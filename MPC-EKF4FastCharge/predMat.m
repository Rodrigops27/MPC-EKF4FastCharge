function [Phi, G, aug] = predMat(A, B, C, D, Np, Nc)
% PREDMAT  Build prediction matrices for y = Phi*x_aug + G*ΔU
% Inputs
%   A (nx×nx), B (nx×m), C (q×nx), D (q×m)  : discrete-time model
%   Np, Nc                                   : prediction/control horizons
% Output
%   Phi (q*Np × (nx+m))                      : stacked output map
%   G   (q*Np × (m*Nc))                      : block Toeplitz from impulse response
%   aug : struct with fields A,B,C and sizes (nx_aug, m, q)

    % ---- sizes & clamps ----
    nx = size(A,1);
    m  = size(B,2);
    q  = size(C,1);
    Np = max(1, round(Np));
    Nc = max(1, min(round(Nc), Np));

    % ---- Δu augmentation ----
    % x_aug = [x; u],  Δu is input
    Abar = [A, B; zeros(m,nx), eye(m)];
    Bbar = [zeros(nx,m); eye(m)];
    Cbar = [C, D];             % feed-through carried in augmented C

    nx_aug = nx + m;

    % ---- preallocate ----
    Phi = zeros(q*Np, nx_aug);
    G   = zeros(q*Np, m*Nc);

    % ---- precompute impulse sequence H_k = Cbar*Abar^k*Bbar, k=0..Np-1 ----
    H = cell(Np,1);
    X = Bbar;                   % Abar^0*Bbar
    for k = 1:Np
        H{k} = Cbar * X;        % q×m
        X    = Abar * X;        % advance: X := Abar^(k)*Bbar
    end

    % ---- build Phi and G without powers or toeplitz ----
    Ap = eye(nx_aug);           % Abar^0
    for i = 1:Np
        Ap      = Abar * Ap;                        % Abar^i
        r       = (i-1)*q + (1:q);
        Phi(r,:) = Cbar * Ap;                       % q×nx_aug

        % place H_{i-j} blocks in row i
        for j = 1:min(i, Nc)
            c = (j-1)*m + (1:m);
            G(r, c) = H{i-j+1};                     % q×m
        end
    end

    % ---- pack aug info ----
    aug = struct('A',Abar,'B',Bbar,'C',Cbar,'nx_aug',nx_aug,'m',m,'q',q);
end
