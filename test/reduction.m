reduction_data;
time = tic;
[reduced_net, fctable, blocked] = QFCA(model, true, 'gurobi');
toc(time);
save('blocked.mat', 'blocked');
save('fctable.mat', 'fctable');
save('reduced_net.mat', 'reduced_net');