function [S, rev, fctable, blocked] = QFCA(S, rev, reduction, varargin)
% QFCA computes the table of flux coupling relations and the list of blocked 
% reactions for a metabolic network specified by its stoichiometric matrix 
% and irreversible reactions and also returns the reduced metabolic network.
%% Usage
% [S_reduced, rev_reduced, fctable, blocked] = QFCA(S, rev, reduction [, solver])
%   - S: the associated sparse stoichiometric matrix
%   - rev: the 0-1 vector with 1's corresponding to the reversible reactions
%   - reduction: logical indicating whether DCE-induced reductions should be
%   carried out or not
%   - solver: the LP solver to be used; the currently available options are
%   either 'linprog' or 'gurobi' with the default value of 'linprog'
%   
%   - S_reduced: the reduced sparse stoichiometric matrix
%   - rev_reduced: the reduced reversibility vector
%   - fctable: the resulting flux coupling matrix
%           * For the choice of entries, we use the F2C2 convention for the
%           sake of compatibility. The meaning of the entry (i, j) is:
%           0 - uncoupled reactions
%           1 - fully coupled reactions
%           2 - partially coupled reactions
%           3 - reaction i is directionally coupled to reaction j
%           4 - reaction j is directionally coupled to reaction i
%   - blocked: the 0-1 vector with 1's corresponding to the blocked reactions
    %% setting up the LP solver
    if ~isempty(varargin)
        solver = varargin{1};
    else
        solver = 'linprog';
    end
    %% identifying the blocked reactions and removing them from the network
    t1 = cputime;
    S = unique(S, 'rows', 'stable');
    [S, rev, blocked] = blockedReac(S, rev, solver);
    % aggregating all the isozymes
    [~, reacNum, duplicates] = unique([S.', rev], 'rows', 'stable');
    duplicates = duplicates.';
    S = S(:, reacNum);
    rev = rev(reacNum);
    fullCouplings = reacNum(duplicates);
    % removing the newly blocked reactions
    [S, rev, newlyBlocked] = blockedReac(S, rev, solver);
    reacNum(newlyBlocked == 1) = [];
    t2 = cputime;
    fprintf('Identifying the blocked reactions and removing them from the network: %.3f\n', t2-t1);
    t1 = t2;
    [m, n] = size(S);
    fprintf('Reduced number of:\n\tmetabolites = %d;\treactions = %d;\tnonzero elements = %d\n', ...
        m, n, nnz(S));
    %% identifying the fully coupled pairs of reactions
    % finding the trivial full coupling relations
    flag = true;
    while flag
        flag = false;
        for i = m:-1:1
            if i <= size(S, 1)
                nzcols = find(S(i, :));
                % check to see if the i-th row of S has only 2 nonzero elements
                if length(nzcols) == 2
                    [S, rev] = mergeFullyCoupled(S, rev, nzcols(1), nzcols(2), ...
                        -S(i, nzcols(1))/S(i, nzcols(2)));
                    fullCouplings(fullCouplings == reacNum(nzcols(2))) = ...
                        reacNum(nzcols(1));
                    reacNum(nzcols(2)) = [];
                    flag = true;
                end
            end
        end
    end
    % finding the rest of full coupling relations
    n = size(S, 2);
    [Q, R, ~] = qr(S.');
    tol = norm(S, 'fro')*eps(class(S));
    % Z is the kernel of the stoichiometric matrix
    Z = Q(:, sum(abs(diag(R)) > tol)+1:n);
    X = tril(Z*Z.');
    Y = diag(diag(X).^(-1/2));
    X = Y*X*Y;
    for i = n:-1:2
        % j is the candidate reaction to be fully coupled to reaction i
        [M, j] = max(abs(X(i, 1:i-1)));
        % this is in fact cauchy-schwarz inequality
        if M > 1 - tol
            % c is the full coupling coefficient
            %c = sign(X(i, j))*Y(j, j)/Y(i, i);
            %if norm(Z(i, :) - c*Z(j, :)) < tol
            [S, rev] = mergeFullyCoupled(S, rev, j, i, sign(X(i, j))*Y(j, j)/Y(i, i));
            fullCouplings(fullCouplings == reacNum(i)) = reacNum(j);
            reacNum(i) = [];
            %end
        end
    end
    S(:, rev == -1) = -S(:, rev == -1);
    rev(rev == -1) = 0;
    t2 = cputime;
    fprintf('Finding the full coupling relations: %.3f\n', t2-t1);
    t1 = t2;
    [m, n] = size(S);
    fprintf('Reduced number of:\n\tmetabolites = %d;\treactions = %d;\tnonzero elements = %d\n', ...
        m, n, nnz(S));
    %% computing the set of fully reversible reactions 
    [~, ~, prev] = blockedReac(S(:, rev == 1), rev(rev == 1), solver);
    % marking the Frev set by 2 in the rev vector
    rev(rev == 1) = 2 - prev;
    t2 = cputime;
    fprintf('Correcting the reversibility types: %.3f\n', t2-t1);
    t1 = t2;
    %% QFCA finds the coupling coefficients
    k = n;
    reacs = 1:n;
    reactions = false(n, 1);
    A = zeros(n);
    for i = n:-1:1
        if rev(i) ~= 2
            %% Irev ---> Irev couplings
            result = directionallyCoupled(S, rev, i, solver);
            dcouplings = result.x(m+1:end) < -0.5;
            dcouplings(i) = false;
            if any(dcouplings)
                A(reacs(rev == 0), i) = 3*dcouplings(rev == 0);
                A(i, reacs(rev == 0)) = 4*dcouplings(rev == 0).';
                dcouplings(i) = true;
                % correcting the reversibility conditions
                if rev(i) == 1
                    rev(i) = 0;
                    if ~reduction && result.x(m+i) < 0
                        S(:, i) = -S(:, i);
                    end
                end
                % inferring by the transitivity of directional coupling relations
                A(reacs(rev == 1), i) = max(A(reacs(rev == 1), ...
                    reacs(dcouplings)), [], 2);
                A(i, reacs(rev == 1)) = max(A(reacs(dcouplings), ...
                    reacs(rev == 1)), [], 1);
                %% Prev ---> Irev couplings
                if any(A(reacs(rev == 1), i) == 0)
                    coupled = false(n, 1);
                    [Q, R, ~] = qr(transpose(S(:, ~dcouplings)));
                    tol = norm(S(:, ~dcouplings), 'fro')*eps(class(S));
                    Z = Q(:, sum(abs(diag(R)) > tol)+1:end);
                    coupled(~dcouplings & rev == 1 & A(reacs, i) == 0) = ...
                        vecnorm(Z(rev(~dcouplings) == 1 & A(reacs(~dcouplings), ...
                        i) == 0, :), 2, 2) < tol;
                    A(reacs(coupled), i) = 3;
                    A(i, reacs(coupled)) = 4;
                    % -1 indicates an uncoupled pair for remembering to skip 
                    % it without any need for further double check later
                    A(reacs(~coupled & rev == 1 & A(reacs, i) == 0), ...
                        reacs(dcouplings)) = -1;
                end
                % metabolic network reduction
                if reduction
                    c = S.'*result.x(1:m);
                    S = S + repmat(S(:, i), 1, n)*spdiags(-c/c(i), 0, n, n);
                    [S, rev] = mergeFullyCoupled(S, rev, i, i, 1);
                    reacs(i) = [];
                    [m, n] = size(S);
                end
                reactions(i) = true;
                for j = i+1:k
                    if reactions(j)
                        if all(A(i, reacs(rev == 0 & ~reactions(reacs))) == ...
                                A(j, reacs(rev == 0 & ~reactions(reacs))))
                            A(i, j) = 2;
                            A(j, i) = 2;
                        elseif all(A(i, reacs(rev == 0 & ~reactions(reacs))) ...
                                <= A(j, reacs(rev == 0 & ~reactions(reacs))))
                            A(i, j) = 3;
                            A(j, i) = 4;
                        elseif all(A(i, reacs(rev == 0 & ~reactions(reacs))) ...
                                >= A(j, reacs(rev == 0 & ~reactions(reacs))))
                            A(i, j) = 4;
                            A(j, i) = 3;
                        end
                    end
                end
            end
        end
    end
    % the usage of -1 was temporary and we return to our earlier convention
    A(A == -1) = 0;
    A(logical(eye(k))) = 1;
    t2 = cputime;
    fprintf('Finding the directional and partial coupling relations: %.3f\n', t2-t1);
    t1 = t2;
    %% postprocessing to fill in the flux coupling table for the original 
    % metabolic network from the flux coupling relations for the reduced one
    map = repmat(fullCouplings.', k, 1) == repmat(reacNum, 1, length(duplicates));
    fctable = map.'*A*map;
    t2 = cputime;
    fprintf('Inferring by the transitivity of full coupling relations: %.3f\n', t2-t1);
    t1 = t2;
    %% reaction pairs that become blocked after merging isozymes are fully coupled
    for i = 1:duplicates(end)
        blockedAfterMerging = find(duplicates == i);
        if (length(blockedAfterMerging) > 2 || ...
                (length(blockedAfterMerging) == 2 && newlyBlocked(i) == 0))
            fctable(duplicates == i, :) = 0;
            fctable(:, duplicates == i) = 0;
        elseif length(blockedAfterMerging) == 2 && newlyBlocked(i) == 1
            fctable(blockedAfterMerging(1), blockedAfterMerging(2)) = 1;
            fctable(blockedAfterMerging(2), blockedAfterMerging(1)) = 1;
        end
    end
    fctable(logical(eye(size(fctable)))) = 1;
    fprintf('Metabolic network reductions postprocessing: %.3f\n', cputime-t1);
    fprintf('Reduced number of:\n\tmetabolites = %d;\treactions = %d;\tnonzero elements = %d\n', ...
        m, n, nnz(S));
end