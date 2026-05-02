n = 75;
m = 4;
p = 0.5;

opts = struct();
opts.seed = 0;
opts.verbose = 0;
opts.maxit = 10;
opts.max_node = 100;


for k = 1:m
    filename = "Random_graph\Erdos_n" + num2str(n) + "_p" + num2str(p) + "_test" + num2str(k);
    load(filename + '.mat', 'A');
    t_start = tic;
    [t, out] = DiSOS_clique_BnB(A, opts);
    t_end = toc(t_start);
    
    disp('--------------------------------')
    save(filename + '_bnb.mat', 'out');
    disp(['time = ', num2str(t_end), ' s, iter = ', num2str(out.iter)])
    disp(['LB = ', num2str(-out.lb_vec(end)), ', UB = ', num2str(-out.ub_vec(end))]) 
end

for k = 1:m
    filename = "Random_graph\Erdos_n" + num2str(n) + "_p" + num2str(p) + "_test" + num2str(k);
    load(filename + '.mat', 'A');
    disp('--------------------------------')
    benchmark(A);   
end


function benchmark(A)
n = size(A, 1);
I = eye(n);
J = ones(n);
Ac = J - I - A;

cvx_begin quiet
cvx_solver MOSEK
    variable t 
    variable N(n, n) symmetric
    minimize(t)
    t*(I+Ac) - J - N == semidefinite(n);
    N >= 0;
cvx_end
disp(['Schrijver theta number (SDP): ', num2str(cvx_optval)])

% integer programming
[I, J] = find(triu(Ac, 1));
cvx_begin quiet
cvx_solver MOSEK
    variable x(n, 1) binary
    maximize(sum(x))
    subject to
        x(I) + x(J) <= 1;
cvx_end
disp(['Clique number (Integer progarmming): ', num2str(cvx_optval)])
end