function [S, rev, fctable, blocked] = QFCA(S, rev, solver, varargin)
% QFCA computes the table of the flux coupling relations and the list of 
% blocked reactions for a metabolic network specified by its stoichiometric 
% matrix and irreversible reactions and also returns the reduced network.
%
% [S_reduced, rev_reduced, fctable, blocked] = QFCA(S, rev, solver [, tol])
%   - S: the sparse stoichiometric matrix
%   - rev: the 0-1 vector with 1's corresponding to reversible reactions
%   - solver: the LP solver to be used; the currently available options are
%   'linprog' and 'gurobi'
%   - tol: the tolerance level with the default value 10e-6
%   - S_reduced: the reduced sparse stoichiometric matrix
%   - rev_reduced: the reduced reversibility vector
%   - fctable: the resulting flux coupling matrix
%           * For the choice of entries we use the F2C2 convention for the
%           sake of compatibility. The meaning of entry (i, j) is:
%           0 - uncoupled reactions
%           1 - fully coupled reactions
%           2 - partially coupled reactions
%           3 - reaction i is directionally coupled to j
%           4 - reaction j is directionally coupled to i
%   - blocked: a 0-1 vector with 1's corresponding to blocked reactions
    %% setting the zero tolerance parameter
    if ~isempty(varargin)
        tol = varargin{1};
    else
        tol = 1e-6;
    end
    %% identifying the blocked reactions and removing them from the network
    t1 = cputime;
    [S, rev, blocked] = blockedReac(S, rev, solver, tol);
    % aggregating all isozymes
    S = unique(S, 'rows');
    [~, reacNum, duplicates] = unique([S.', rev], 'rows', 'stable');
    duplicates = duplicates.';
    S = S(:, reacNum);
    rev = rev(reacNum);
    fullCouplings = reacNum(duplicates);
    % removing the newly blocked reactions
    [S, rev, newBlocked] = blockedReac(S, rev, solver, tol);
    reacNum(newBlocked == 1) = [];
    t2 = cputime;
    fprintf('Identifying the blocked reactions and removing them from the network: %.3f\n', t2-t1);
    t1 = t2;
    [m, n] = size(S);
    fprintf('\nReduced number of metabolites = %d;\t', m);
    fprintf('Reduced number of reactions = %d\n', n);
    %% identifying the fully coupled pairs of reactions
    % finding the trivial full couplings
    flag = true;
    while flag
        flag = false;
        for i = m:-1:1
            if i <= size(S, 1)
                nzcols = find(S(i, :));
                % check if the i-th row of S has only 2 nonzero elements
                if length(nzcols) == 2
                    [S, rev] = mergeFullyCoupled(S, rev, nzcols(1), nzcols(2), -S(i, nzcols(1))/S(i, nzcols(2)));
                    fullCouplings(fullCouplings == reacNum(nzcols(2))) = reacNum(nzcols(1));
                    reacNum(nzcols(2)) = [];
                    flag = true;
                end
            end
        end
    end
    % finding the rest of full couplings
    n = size(S, 2);
    [Q, R, ~] = qr(S.');
    s = abs(diag(R));
    r = sum(s > tol);
    Z = Q(:, r+1:n);
    X = Z*Z.';
    for i = n:-1:2
        cand = find(X(i, i)*diag(X(1:i-1, 1:i-1))./X(1:i-1, i).^2 < 1 + tol);
        % cand is the list of candidate reactions to be fully coupled to i
        if ~isempty(cand)
            for j = cand.'
                c = sqrt(X(i, i)/X(j, j));
                % c is the absolute value of the full coupling coefficient
                if norm(Z(i, :) - c*Z(j, :)) < tol
                    [S, rev] = mergeFullyCoupled(S, rev, j, i, c);
                    fullCouplings(fullCouplings == reacNum(i)) = reacNum(j);
                    reacNum(i) = [];
                    break;
                elseif norm(Z(i, :) + c*Z(j, :)) < tol
                    [S, rev] = mergeFullyCoupled(S, rev, j, i, -c);
                    fullCouplings(fullCouplings == reacNum(i)) = reacNum(j);
                    reacNum(i) = [];
                    break;
                end
            end
        end
    end
    S(:, rev == -1) = -S(:, rev == -1);
    rev(rev == -1) = 0;
    t2 = cputime;
    fprintf('Finding the full coupling relations: %.3f\n', t2-t1);
    t1 = t2;
    [m, n] = size(S);
    fprintf('\nReduced number of metabolites = %d;\t', m);
    fprintf('Reduced number of reactions = %d\n', n);
    %% correcting the reversibility conditions 
    % computing the set of fully reversible reactions 
    [~, ~, prev] = blockedReac(S(:, rev == 1), rev(rev == 1), solver, tol);
    % marking the Frev set by 2 in the rev vector
    rev(rev == 1) = 2 - prev;
    t2 = cputime;
    fprintf('Correcting the reversibility types: %.3f\n', t2-t1);
    t1 = t2;
    %% QFCA finds the coupling coefficients
    reacs = 1:n;
    A = zeros(n);
    for i = n:-1:1
        if rev(i) ~= 2
            %% Irev ---> Irev couplings
            result = directionallyCoupled(S, rev, i, solver);
            dcouplings = result.x(m+1:end) < -0.5;
            dcouplings(i) = false;
            A(reacs, reacs(i)) = 3*dcouplings;
            A(reacs(i), reacs) = 4*dcouplings.';
            %% Prev ---> Irev couplings
            if any(dcouplings)
                dcouplings(i) = true;
                % inferring by the transitivity of directional couplings
                A(reacs(rev == 1), reacs(i)) = max(A(reacs(rev == 1), reacs(dcouplings)), [], 2);
                A(reacs(i), reacs(rev == 1)) = max(A(reacs(dcouplings), reacs(rev == 1)), [], 1);
                couplings = ~dcouplings;
                B = speye(sum(couplings));
                B = B(:, rev(couplings) == 1 & A(reacs(couplings), reacs(i)) == 0);
                % the columns of B correspond to the candidate reactions
                if ~isempty(B)
                    X = mldivide(transpose(S(:, couplings)), B);
                    coupled = false(n, 1);
                    coupled(couplings & rev == 1 & A(reacs, reacs(i)) == 0) = all(abs(transpose(S(:, couplings))*X-B) < tol, 1);
                    A(reacs(coupled), reacs(i)) = 3;
                    A(reacs(i), reacs(coupled)) = 4;
                    % -1 indicates an uncoupled pair for remembering to skip 
                    % it without any need for further double check later
                    A(reacs(~coupled & rev == 1 & A(reacs, reacs(i)) == 0), reacs(dcouplings)) = -1;
                end
                % metabolic network reduction
                c = S.'*result.x(1:m);
                S = S + repmat(S(:, i), 1, n)*spdiags(-c/c(i), 0, n, n);
                [S, rev] = mergeFullyCoupled(S, rev, i, i, 1);
                reacs(i) = [];
                [m, n] = size(S);
            end
        end
    end
    % the usage of -1 was temporary and we return to our earlier convention
    A(A == -1) = 0;
    k = length(A);
    A(logical(eye(k))) = 1;
    t2 = cputime;
    fprintf('Finding the directional and partial couplings: %.3f\n', t2-t1);
    t1 = t2;
    %% postprocessing to fill in the flux coupling table for the original 
    % network from the flux coupling relations for the reduced one
    for i = 1:k-1
        if all(reacs ~= i)
            for j = i+1:k
                if all(reacs ~= j)
                    if all(A(i, reacs) == A(j, reacs))
                        A(i, j) = 2;
                        A(j, i) = 2;
                    elseif all(A(i, reacs) <= A(j, reacs))
                        A(i, j) = 3;
                        A(j, i) = 4;
                    elseif all(A(i, reacs) >= A(j, reacs))
                        A(i, j) = 4;
                        A(j, i) = 3;
                    end
                end
            end
        end
    end
    t2 = cputime;
    fprintf('Metabolic network reductions postprocessing: %.3f\n', t2-t1);
    t1 = t2;
    %% inferring by the transitivity of full couplings
    l = length(duplicates);
    map = repmat(fullCouplings.', k, 1) == repmat(reacNum, 1, l);
    fctable = map.'*A*map;
    fctable(logical(eye(l))) = 0;
    for i = unique(duplicates)
        blockedAfterMerging = find(duplicates == i);
        % pairs that become blocked after merging isozymes are fully coupled
        if length(blockedAfterMerging) > 2 || (length(blockedAfterMerging) == 2 && newBlocked(i) == 0)
            fctable(duplicates == i, :) = 0;
            fctable(:, duplicates == i) = 0;
        elseif length(blockedAfterMerging) == 2 && newBlocked(i) == 1
            fctable(blockedAfterMerging(1), blockedAfterMerging(2)) = 1;
            fctable(blockedAfterMerging(2), blockedAfterMerging(1)) = 1;
        end
    end
    fctable = fctable + eye(l);
    fprintf('Inferring by the transitivity of full couplings: %.3f\n', cputime-t1);
    fprintf('\nReduced number of metabolites = %d;\t', m);
    fprintf('Reduced number of reactions = %d\n', n);
end