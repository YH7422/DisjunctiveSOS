n = 3;
x = sdpvar(n, 1);
opts = struct();
opts.init = 1;
opts.method = 0;
opts.verbose = 0;
opts.eps = 1e-6;
opts.maxit = 300;
opts.dq = 2;

for m = 1:6
    d = 2*(2*m+1);
    p = x(1)^(2*m+1) * x(3)^(2*m+1) + (x(2)^2 * x(3)^(2*m-1) - x(1)^(2*m+1) - x(1)*x(3)^(2*m))^2;
    t_start = tic;
    [t, out] = DiSOS_BnB_SD(x, p, n, d, opts);
    t_end = toc(t_start);
    disp(['Stengle-', num2str(m), ': time = ', num2str(t_end), ' s, iter = ', num2str(out.iter)])
    save(['Stengle', num2str(m), '_bnb_sd_init', num2str(opts.init)], 'out');
    % save the sequence of lower and upper bounds
end