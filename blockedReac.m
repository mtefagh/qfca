function [S, rev, blocked] = blockedReac(S, rev, solver)
% blockedReac finds the blocked reactions and removes them from the network
    [m, n] = size(S);
    blocked = zeros(n, 1);
    % identifying the blocked irreversible reactions
    result = directionallyCoupled(S, rev, 0, solver);
    blocked(result.x(m+1:end) < -0.5) = 1;
    % identifying the blocked reversible reactions
    [Q, R, ~] = qr(transpose(S(:, blocked == 0)));
    % setting up the zero-tolerance parameter
    tol = norm(S(:, blocked == 0), 'fro')*eps(class(S));
    Z = Q(:, sum(abs(diag(R)) > tol)+1:end);
    blocked(blocked == 0 & rev == 1) = vecnorm(Z(rev(blocked == 0) == 1, :), ...
        2, 2) < tol;
    % removing the blocked reactions from the metabolic network
    S(:, blocked == 1) = [];
    rev(blocked == 1) = [];
end