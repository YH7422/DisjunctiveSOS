m = 4;
Q_array = cell(m, 1);

Q_array{1} = [1 0 1 1 0;
0 1 0 1 1;
1 0 1 0 1;
1 1 0 1 0;
0 1 1 0 1];

Q_array{3} = [-14 -15 -16 0 0;
-15 -14 -12.5 -22.5 -15;
-16 -12.5 -10 -26.5 -16;
0 -22.5 -26.5 0 0;
0 -15 -16 0 -14];

Q_array{4} = [0.9044 0.1054 0.5140 0.3322 0;
0.1054 0.8715 0.7385 0.5866 0.9751;
0.5140 0.7385 0.6936 0.5368 0.8086;
0.3322 0.5866 0.5368 0.5633 0.7478;
0 0.9751 0.8086 0.7478 1.2932];

Q_array{2} = [1 0 0 0 0 0 1 1 1 1 1 1;
0 1 0 0 1 1 0 0 1 1 1 1;
0 0 1 1 0 1 0 1 0 1 1 1;
0 0 1 1 1 0 1 0 1 0 1 1;
0 1 0 1 1 0 1 1 0 1 0 1;
0 1 1 0 0 1 1 1 1 0 0 1;
1 0 0 1 1 1 1 0 0 1 1 0;
1 0 1 0 1 1 0 1 1 0 1 0;
1 1 0 1 0 1 0 1 1 1 0 0;
1 1 1 0 1 0 1 0 1 1 0 0;
1 1 1 1 0 0 1 1 0 0 1 0;
1 1 1 1 1 1 0 0 0 0 0 1];


opts = struct();
opts.seed = 0;
opts.verbose = 0;
opts.max_node = 1000;

for i = 1:4
    Q = Q_array{i};
    n = size(Q, 1);
    t_start = tic;
    [t, out] = DiSOS_copositive_BnB(Q, opts);
    t_end = toc(t_start);
    disp(['Test ', num2str(i), ': time = ', num2str(t_end), ' s, iter = ', num2str(out.iter)])
    save(['copositive', num2str(i), '_bnb'], 'out');
    % save the sequence of lower and upper bounds
 
%     figure
%     plot(out.lb_vec, 'Displayname', 'lower bound')
%     hold on 
%     plot(out.ub_vec, 'Displayname', 'upper bound')
%     xlabel('Iterations')
%     ylabel('Objective value')
%     legend('Location', 'best')
end