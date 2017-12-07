function [S, rev] = mergeFullyCoupled(S, rev, i, j, c)
    %% merging the fully coupled reactions
    S(:, i) = S(:, i) + c*S(:, j);
    S(:, j) = [];
    % deleting zero rows from the stoichiometric matrix
    S = S(any(logical(S), 2), :);
    % deleting the reaction from the rev vector
    if rev(j) ~= 1
        if c > 0
            rev(i) = rev(j);
        else
            rev(i) = -1 - rev(j);
        end
    end
    rev(j) = [];
end