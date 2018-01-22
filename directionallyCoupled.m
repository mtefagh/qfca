function result = directionallyCoupled(S, rev, i, solver)
% directionallyCoupled finds all the directionally coupled reactions to i
% the currently available options for the LP solver are linprog and gurobi
    [m, n] = size(S);
    irevIndex = [m+1:m+i-1, m+i+1:m+n];
    irevIndex = irevIndex(rev([1:i-1, i+1:n]) == 0);
    model.obj = zeros(m+n, 1);
    model.obj(irevIndex) = 1;
    model.A = [S', -speye(n)];
    model.sense = repmat('=', n, 1);
    model.sense(rev == 0) = '<';
    model.rhs = zeros(n, 1);
    model.lb = [-Inf(m, 1); zeros(n, 1)];
    model.lb(m+i) = -Inf;
    model.lb(irevIndex) = -1;
    model.ub = [Inf(m, 1); zeros(n, 1)];
    model.ub(m+i) = Inf;
    if strcmp(solver, 'gurobi')
        params.outputflag = 0;
        result = gurobi(model, params);
        if ~strcmp(result.status, 'OPTIMAL')
            warning('Optimization is unstable!');
            disp(result.status);
        end
    elseif strcmp(solver, 'linprog')
        problem.f = model.obj;
        problem.Aineq = model.A(rev == 0, :);
        problem.bineq = model.rhs(rev == 0);
        problem.Aeq = model.A(rev ~= 0, :);
        problem.beq = model.rhs(rev ~= 0);
        problem.lb = model.lb;
        problem.ub = model.ub;
        problem.solver = 'linprog';
        problem.options = optimset('Display', 'off');
        [result.x, result.objval, result.status, ~] = linprog(problem);
        if result.status ~= 1
            warning('Optimization is unstable!');
            disp(result.status);
        end
    else
        warning('The entered solver is not supported in this version!');
    end
end