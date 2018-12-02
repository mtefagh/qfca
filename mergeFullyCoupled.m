function [S, rev] = mergeFullyCoupled(S, rev, i, j, c)
% mergeFullyCoupled merges the fully coupled pair of reactions (i, j)
    S(:, i) = S(:, i) + c*S(:, j);
    S(:, j) = [];
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