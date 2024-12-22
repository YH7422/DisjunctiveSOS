qvec = [17, 29, 37, 41];

opts = struct();
opts.seed = 0;
opts.verbose = 0;

for q = qvec
    filename = ['Paley_graph\Paley_q', num2str(q), '.csv'];
    A = readmatrix(filename);
    t_start = tic;
    [t, out] = DiSOS_clique_BnB(A, opts);
    t_end = toc(t_start);
    
    disp('--------------------------------')
    disp(['Paley q = ', num2str(q), ': time = ', num2str(t_end), ' s, iter = ', num2str(out.iter)])
    benchmark(A);    
end


function benchmark(A)
q = size(A, 1);

ub = (sqrt(2*q - 1) + 1)/2;
disp(['Theoretical upper bound: ', num2str(ub), '(', num2str(floor(ub)), ')'])

ub = sqrt(q);
disp(['P+N upper bound: ', num2str(ub), '(', num2str(floor(ub)), ')'])

% integer programming
cvx_begin quiet
cvx_solver MOSEK
    variable x(q, 1) binary;
    maximize(sum(x))
    for i = 1 : q-1
        for j = i+1 : q
            if A(i,j) == 0
                x(i) + x(j) <= 1;
            end
        end
    end
cvx_end
disp(['Clique number (Integer progarmming): ', num2str(cvx_optval)])
end