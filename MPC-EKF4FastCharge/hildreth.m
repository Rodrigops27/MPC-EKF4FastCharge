function [DU, lambda, nexec] = hildreth(E, F, M, gamma, lambda0, maxIter)
% Solve:  minimize (1/2) DU' * E * DU + F' * DU
%         s.t.      M * DU <= gamma
%
% Inputs:
%   E       : Hessian (positive definite)
%   F       : gradient
%   M,gamma : inequality constraints
%   lambda0 : initial guess for multipliers
%   maxIter : maximum iterations
%
% Outputs:
%   DU      : optimal decision vector
%   lambda  : optimal multipliers
%   nexec   : number of iterations executed

    nConstr = size(M,1);
    if nargin < 5 || isempty(lambda0)
        lambda = zeros(nConstr,1);
    else
        lambda = lambda0;
    end
    if nargin < 6
        maxIter = 100;
    end
    
    % Precompute H and K (dual system)
    H = M*(E\M');       % nConstr x nConstr
    K = M*(E\F) + gamma; 
    
    % Iterative updates
    for k = 1:maxIter
        lambda_old = lambda;
        for i = 1:nConstr
            w = -(K(i) + H(i,:)*lambda - H(i,i)*lambda(i)) / H(i,i);
            lambda(i) = max(0, w);
        end
        % Check convergence
        if norm(lambda - lambda_old, Inf) < 1e-6
            break;
        end
    end
    nexec = k;
    
    % Recover DU from primal relation
    DU = -E\(F + M'*lambda);
end
