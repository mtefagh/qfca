function [fctable, blocked] = QFCA(S, rev, solver, varargin)
    % setting the zero tolerance parameter
    if ~isempty(varargin)
        tol = varargin{1};
    else
        tol = 1e-6;
    end
    t1 = cputime;
    % identifying the blocked reactions and removing them from the network
    [S, rev, blocked] = blockedReac(S, rev, solver);
    % aggregating all isozymes
    n = size(S, 2);
    duplicates = 1:n;
    reacNum = 1:n;
    for i = n:-1:2
        for j = i-1:-1:1
            if all(S(:, i) == S(:, j)) && rev(i) == rev(j)
                [S, rev] = mergeFullyCoupled(S, rev, i, i, 1);
                duplicates(duplicates == i) = j;
                reacNum(i) = [];
                break;
            end
        end
    end
    % removing the newly blocked reactions
    [S, rev, newBlocked] = blockedReac(S, rev, solver);
    blockedAfterMerging = reacNum(newBlocked == 1);
    reacNum(newBlocked == 1) = [];
    t2 = cputime;
    fprintf('Identifying the blocked reactions and removing them from the network: %f\n', t2-t1);
    t1 = t2;
    % finding the trivial full couplings
    fullCouplings = duplicates;
    flag = true;
    while flag
        flag = false;
        for i = size(S, 1):-1:1
            if i <= size(S, 1)
                nzcols = find(S(i, :));
                if length(nzcols) == 2
                    [S, rev] = mergeFullyCoupled(S, rev, nzcols(1), nzcols(2), -S(i, nzcols(1))/S(i, nzcols(2)));
                    fullCouplings(fullCouplings == reacNum(nzcols(2))) = reacNum(nzcols(1));
                    reacNum(nzcols(2)) = [];
                    flag = true;
                end
            end
        end
    end
    % finding the rest of fully couplings
    [Q, R, ~] = qr(S.');
    s = abs(diag(R));
    r = sum(s > 1e-6);
    Z = Q(:, r+1:size(S, 2));
    X = Z*Z.';
    for i = size(S, 2):-1:2
        cand = find(X(i, i)*diag(X(1:i-1, 1:i-1))./X(1:i-1, i).^2 < 1 + tol);
        if ~isempty(cand)
            for j = cand.'
                c = sqrt(X(i, i)/X(j, j));
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
    fprintf('Finding the full couplings: %f\n', t2-t1);
    t1 = t2;
    % computing the set of fully reversible reactions and mark the Frev set
    % by 2 in rev vector
    [~, ~, prev] = blockedReac(S(:, rev == 1), rev(rev == 1), solver);
    rev(rev == 1) = 2 - prev;
    % correcting the reversibility conditions and mark the Irev set by 0 in
    % rev vector
    [m, n] = size(S);
    for i = 1:n
        if rev(i) == 1
            result = directionallyCoupled(S, rev, i, solver);
            if result.objval < -0.5
                rev(i) = 0;
                if result.x(m+i) < 0
                    S(:, i) = -S(:, i);
                end
            end
        end
    end
    t2 = cputime;
    fprintf('Correcting the reversibility types: %f\n', t2-t1);
    t1 = t2;
    % QFCA finds the coupling coefficients
    reacs = 1:n;
    A = eye(n);
    for i = n:-1:1
        if rev(i) == 0
            % Irev ---> Irev couplings
            result = directionallyCoupled(S, rev, i, solver);
            dcouplings = result.x(m+1:end) < -0.5;
            A(reacs, reacs(i)) = A(reacs, reacs(i)) + 3*dcouplings;
            A(reacs(i), reacs) = A(reacs(i), reacs) + 4*dcouplings.';
            % Prev ---> Irev couplings
            if any(dcouplings)
                dcouplings(i) = true;
                A(reacs(rev == 1), reacs(i)) = max(A(reacs(rev == 1), reacs(dcouplings)), [], 2);
                A(reacs(i), reacs(rev == 1)) = max(A(reacs(dcouplings), reacs(rev == 1)), [], 1);
                couplings = ~dcouplings;
                B = speye(sum(couplings));
                B = B(:, rev(couplings) == 1 & A(reacs(couplings), reacs(i)) == 0);
                if ~isempty(B)
                    X = mldivide(transpose(S(:, couplings)), B);
                    coupled = false(n, 1);
                    coupled(couplings & rev == 1 & A(reacs, reacs(i)) == 0) = all(abs(transpose(S(:, couplings))*X-B) < tol, 1);
                    A(reacs(coupled), reacs(i)) = 3;
                    A(reacs(i), reacs(coupled)) = 4;
                    A(reacs(~coupled & rev == 1 & A(reacs, reacs(i)) == 0), reacs(dcouplings)) = -1;
                end
                % network reduction
                c = S.'*result.x(1:m);
                S = S + repmat(S(:, i), 1, n)*spdiags(-c/c(i), 0, n, n);
                [S, rev] = mergeFullyCoupled(S, rev, i, i, 1);
                reacs(i) = [];
                [m, n] = size(S);
            end
        end
    end
    A(A == 7) = 2;
    A(A == -1) = 0;
    t2 = cputime;
    fprintf('Finding the directional and partial couplings: %f\n', t2-t1);
    t1 = t2;
    % infering by transitivity of fully couplings
    n = length(A);
    for i = 1:n-1
        if all(reacs ~= i)
            for j = i+1:n
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
    % Post-processings to fill in the the coupling table for the original
    % network from the relations for the reduced one
    fctable = eye(length(duplicates));
    for i = 1:n
        for j = i:n
            fctable(fullCouplings == reacNum(i), fullCouplings == reacNum(j)) = A(i, j);
            fctable(fullCouplings == reacNum(j), fullCouplings == reacNum(i)) = A(j, i);
        end
    end
    fctable = fctable - eye(length(duplicates));
    for i = unique(duplicates)
        if sum(duplicates == i) > 1
            fctable(duplicates == i, :) = 0;
            fctable(:, duplicates == i) = 0;
        end
    end
    for i = blockedAfterMerging
        temp = find(duplicates == i);
        if length(temp) == 2
            fctable(temp(1), temp(2)) = 1;
            fctable(temp(2), temp(1)) = 1;
        end
    end
    fctable = fctable + eye(length(duplicates));
    t2 = cputime;
    fprintf('Post-processings: %f\n', t2-t1);
end