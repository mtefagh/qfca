function [S, rev, blocked] = blockedReac(S, rev, solver, tol)
% blockedReac finds the blocked reactions and removes them from the network
    [m, n] = size(S);
    blocked = zeros(n, 1);
    % identifying the blocked irreversible reactions
    result = directionallyCoupled(S, rev, 0, solver);
    blocked(result.x(m+1:end) < -0.5) = 1;
    % identifying the blocked reversible reactions
    B = speye(sum(blocked == 0));
    B = B(:, rev(blocked == 0) == 1);
    X = mldivide(transpose(S(:, blocked == 0)), B);
    blocked(blocked == 0 & rev == 1) = all(abs(transpose(S(:, blocked == 0))*X-B) < tol, 1);
    % removing the blocked reactions
    S(:, blocked == 1) = [];
    rev(blocked == 1) = [];
end