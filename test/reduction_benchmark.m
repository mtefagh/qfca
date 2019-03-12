reduction_data;
load('blocked.mat');
load('reduced_net.mat');
solver = 'gurobi';
LP = zeros(2, 100);
runtime = zeros(2, 100);
errors = true(2, 100);
unblocked = find(blocked == 0);
for i = 1:100
    disp(i);
    core = randsample(unblocked, 150*i);
    reduced_core = zeros(size(core));
    for k = 1:length(reduced_net.rxns)
        temp = split(reduced_net.rxns{k}, ", ");
        for j = 1:length(core)
            if any(strcmp(model.rxns{core(j)}, temp))
                reduced_core(j) = k;
            end
        end
    end
    reduced_core = unique(reduced_core);
    tic;
    [coreInd, numLP] = swiftcore(reduced_net.S, reduced_net.rev, reduced_core, ones(size(reduced_net.rev)), false, solver);
    runtime(1, i) = toc;
    LP(1, i) = numLP; 
    if all(coreInd(reduced_core))
        if all((1:sum(coreInd)).' == swiftcc(reduced_net.S(:, coreInd), reduced_net.rev(coreInd), solver))
            errors(1, i) = false;
        end
    end
    tic;
    [coreInd, numLP] = swiftcore(model.S, model.rev, core, ones(size(model.rev)), false, solver);
    runtime(2, i) = toc;
    LP(2, i) = numLP;
    if all(coreInd(core))
        if all((1:sum(coreInd)).' == swiftcc(model.S(:, coreInd), model.rev(coreInd), solver))
            errors(2, i) = false;
        end
    end
end
fprintf('SWIFTCORE w/ reduction returns wrong answers in %d%% of cases!\n', round(100*sum(errors(1,:))/size(errors,2)));
fprintf('SWIFTCORE w/o reduction returns wrong answers in %d%% of cases!\n', round(100*sum(errors(2,:))/size(errors,2)));
fprintf('switftcore w/ reduction was computed in %d%% of switftcore w/o reduction runtime.\n', round(mean(runtime(1, :)./runtime(2, :))*100));
% drawing scatterplots for the comparison of the number of LPs
figure();
plot(150*(1:100),LP(1,:),'o',150*(1:100),LP(2,:),'*');
legend({'SWIFTCORE w/ reduction','SWIFTCORE w/o reduction'},'Location','southeast');
xlabel('number of core reactions');
ylabel('number of solved LPs');
savefig('LP');
% drawing scatterplots for the comparison of the runtimes
figure();
plot(150*(1:100),runtime(1,:),'o',150*(1:100),runtime(2,:),'*');
legend({'SWIFTCORE w/ reduction','SWIFTCORE w/o reduction'},'Location','northeast');
xlabel('number of core reactions');
ylabel('runtime in seconds');
savefig('runtime');