qvec = [17, 29, 37, 41, 61];

opts = struct();
opts.seed = 0;
opts.maxit = 10;
opts.verbose = 0;

for q = qvec
    filename = ['Paley_graph\Paley_q', num2str(q), '.csv'];
    A = readmatrix(filename);
    [t, out] = DiSOS_clique_BnB(A, opts);
    
    figure
    plot(-out.ub_vec, 'Displayname', 'lower bound')
    hold on 
    plot(-out.lb_vec, 'Displayname', 'upper bound')
    benchmark(A);    
    xlabel('Iterations')
    legend('Location', 'best')
end


function benchmark(A)
q = size(A, 1);

ub = (sqrt(2*q - 1) + 1)/2;
disp(['Theoretical upper bound: ', num2str(ub), '(', num2str(floor(ub)), ')'])
yline(ub, 'DisplayName', 'Theoretical upper bound', 'Color', [0 0.5 0.5])
hold on

ub = sqrt(q);
disp(['P+N upper bound: ', num2str(ub), '(', num2str(floor(ub)), ')'])
yline(ub, 'DisplayName', 'P+N upper bound', 'Color', [0.5 0 0.5])

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
yline(cvx_optval, 'DisplayName', 'Clique number (Integer progarmming)')
end