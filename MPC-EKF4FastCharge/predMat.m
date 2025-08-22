function [Phi, G] = predMat(A, B, C, D, Np, Nc)
% PREDMAT  Build prediction matrices for y = Phi*x + G*DU
% A,B,C,D : discrete-time, augmented (e.g., x_aug=[x;u])
% Np,Nc   : prediction & control horizons
%
% Works for scalar/multi-output C, and any state size.
% Robust to ill-conditioning: avoids large explicit kron/Toeplitz builds.

nx = size(A,1);
ny = size(C,1);

% Clamp Nc to Np
Nc = min(Nc, Np);

% Preallocate
Phi = zeros(ny*Np, nx);
G   = zeros(ny*Np, Nc);

% Running powers of A (safe single pass)
A_pow = eye(nx);

% Column blocks for G (impulse response of (A,B) mapped by C)
%   y[k+i|k] = C*A^i x[k] + sum_{j=0}^{i-1} C*A^{i-1-j}*B * DU[k+1+j] +
% (C*A^{i-1}*B + D)*u_k  (if not augmented)  TBD
% In Δu + augmented-state design, D is typically 0; we keep it anyway.
gcol = cell(Nc,1);
for j = 1:Nc
    gcol{j} = zeros(ny*Np,1);
end

for i = 1:Np
    % Top block row index for this i
    r = (i-1)*ny + (1:ny);

    % Phi block row
    A_pow = A_pow * A;        % A^i
    Phi(r,:) = C * A_pow;     % C*A^i

    % Fill G block row with C*A^{i-1-j}*B terms
    Ap = eye(nx);             % A^{i-1-j} as we move across j
    for j = 1:min(i,Nc)
        if j==1
            Ap = eye(nx);     % when i-1-j = i-2,... handled below
        end
        % For position (i,j): exponent is (i-1)-(j-1) = i-j
        % We'll compute A^(i-j) by multiplying once per loop
        if j==1
            Aij = A_pow / A;  % A^(i-1)
        else
            Aij = Aij / A;    % step back one power: A^(i-j)
        end

        g = C * (Aij * B);    % ny x 1 (nu=1)
        G(r,j) = g;
    end
end

% If D ≠ 0 and you are NOT using input augmentation, the prediction law
% would include a bias term that depends on u_k. In your Δu+augmented
% design we expect D = 0 (or absorbed in the augmentation), so nothing to add.

useAug = (nargin>=4) && any(abs(D(:))>1e-10);
if useAug
    % Augmented (recommended)
    n = size(A,1); m = size(B,2); q = size(C,1);
    Abar = [A, B; zeros(m,n), eye(m)];
    Bbar = [zeros(n,m); eye(m)];
    Cbar = [C, D];

    % Phi rows: Cbar*Abar^i * [xk; uk]
    Phi = zeros(Np*q, n+m);
    G   = zeros(Np*q, Nc*m);
    Ap  = eye(n+m);
    for i = 1:Np
        Ap  = Abar*Ap;
        Phi((i-1)*q+1:i*q,:) = Cbar*Ap;
        Aj = eye(n+m);
        for j = 1:min(i,Nc)
            G((i-1)*q+1:i*q, (j-1)*m+1:j*m) = Cbar*Aj*Bbar;
            Aj = Abar*Aj;
        end
    end
else
    % Non-augmented (only if D≈0)
    n = size(A,1); m = size(B,2); q = size(C,1);
    Phi = zeros(Np*q, n);
    G   = zeros(Np*q, Nc*m);
    Ap  = eye(n);
    for i = 1:Np
        Ap  = A*Ap;
        Phi((i-1)*q+1:i*q,:) = C*Ap;
        Aj = eye(n);
        for j = 1:min(i,Nc)
            G((i-1)*q+1:i*q, (j-1)*m+1:j*m) = C*Aj*B;
            Aj = A*Aj;
        end
    end
end


end
