name_vec = ["Motzkin"; "Robinson1"; "Robinson2"; "Choi-Lam1"; "Choi-Lam2"; ...
    "Horn"; "Lax"; "Schmudgen"];
n_vec = [3, 3, 4, 4, 3, 5, 5, 3];
d_vec = [6, 6, 4, 4, 6, 4, 4, 6];
n_max = max(n_vec);
x = sdpvar(n_max, 1);
p_vec = [x(1)^4*x(2)^2 + x(1)^2*x(2)^4 - 3*x(1)^2*x(2)^2*x(3)^2 + x(3)^6;
    x(1)^6 + x(2)^6 + x(3)^6 - (x(1)^4*x(2)^2 + x(1)^2*x(2)^4 + x(1)^4*x(3)^2 + x(1)^2*x(3)^4 + ...
    x(2)^4*x(3)^2 + x(2)^2*x(3)^4) + 3*x(1)^2*x(2)^2*x(3)^2;
    x(1)^2*(x(1)-x(4))^2 + x(2)^2*(x(2)-x(4))^2 + x(3)^2*(x(3)-x(4))^2 + 2*x(1)*x(2)*x(3)*(x(1)+x(2)+x(3)-2*x(4));
    x(1)^2*x(2)^2 + x(1)^2*x(3)^2 + x(2)^2*x(3)^2 + x(4)^4 - 4*x(1)*x(2)*x(3)*x(4);
    x(1)^4*x(2)^2 + x(2)^4*x(3)^2 + x(3)^4*x(1)^2 - 3*x(1)^2*x(2)^2*x(3)^2;
    sum(x(1:5).^2)^2 - 4*(x(1)^2*x(2)^2 + x(2)^2*x(3)^2 + x(3)^2*x(4)^2 + x(4)^2*x(5)^2 + x(5)^2*x(1)^2);
    (x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))*(x(1)-x(5)) + (x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))*(x(2)-x(5)) + ...
    (x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))*(x(3)-x(5)) + (x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))*(x(4)-x(5)) + ...
    (x(5)-x(1))*(x(5)-x(2))*(x(5)-x(3))*(x(5)-x(4));
    200*(x(1)^3-4*x(1)*x(3)^2)^2 + 200*(x(2)^3-4*x(2)*x(3)^2)^2 + ...
    (x(2)^2-x(1)^2)*x(1)*(x(1)+2*x(3))*(x(1)^2-2*x(1)*x(3)+2*x(2)^2-8*x(3)^2)];

test_index = (1:length(n_vec));

opts = struct();
opts.init = 0;
opts.method = 0;
opts.verbose = 0;
opts.eps = 1e-6;
opts.maxit = 300;

for i = test_index
    name = name_vec(i);
    n = n_vec(i);
    d = d_vec(i);
    p = p_vec(i);
    t_start = tic;
    [t, out] = DiSOS_BnB_SD(x(1:n), p, n, d, opts);
    t_end = toc(t_start);
    disp(['Test ', num2str(i), ': time = ', num2str(t_end), ' s, iter = ', num2str(out.iter)])
    save([char(name), '_bnb_sd_init', num2str(opts.init)], 'out');
    % save the sequence of lower and upper bounds
end
