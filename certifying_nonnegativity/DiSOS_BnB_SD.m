% A branch-and-bound method aimimg to solve
% min p(x) s.t. x \in S^{n-1}.
% We approximate the SOS cone with
% p(Vx) has nonnegative coefficients on each simplicial cone. 

function [t, out] = DiSOS_BnB_SD(x, p, n, d, opts)
if ~isfield(opts, 'maxit'); opts.maxit = 100; end
if ~isfield(opts, 'eps'); opts.eps = 1e-6; end
% The degree of cuts
if ~isfield(opts, 'init'); opts.init = 0; end
% init = 0, n+1 initial nodes generated from regular simplex
% init = 1, 2^(n-1) initial nodes generated from orthants
if ~isfield(opts, 'method'); opts.method = 0; end
% method = 0, p(V(x.^2)) is sos
% method = 1, p(Vx) has nonnegative coefficients 
if ~isfield(opts, 'solver'); opts.solver = 'MOSEK'; end
if ~isfield(opts, 'verbose'); opts.verbose = 1; end
% verbose = 0 no output
% verbose = 1 basic output
% verbose = 2 detailed output
if ~isfield(opts, 'time'); opts.time = 0; end
if ~isfield(opts, 'seed'); opts.seed = 0; end


% fix the seed to reproduce the result
ss = RandStream('mt19937ar', 'Seed', opts.seed); 
RandStream.setGlobalStream(ss);

if opts.time
    time0 = tic;
    out.tvec = [];
end
if opts.init == 0
    V = [sqrt(1+1/n)*eye(n) - n^(-3/2)*(sqrt(n+1)-1)*ones(n), -ones(n, 1)/sqrt(n)];
    L = -inf;
    ub0 = zeros(n+1, 1);
    for i = 1:n+1
        ub0(i) = replace(p, x, V(:, i));
    end
    U = min(ub0);
    
    % initial nodes
    H = MinHeap_BnB(opts.maxit*2);
    for i = 1:n+1
        Vi = V;
        Vi(:, i) = [];
        lb = get_lower_bound(x, p, n, d, Vi, opts);
        L = max(L, lb);
        center = mean(Vi, 2);
        ubi = ub0;
        ubi(i) = [];
        ub = min(min(ubi), replace(p, x, center));
        U = min(U, ub);
        node = Node([], [], lb, ub, Vi);
        H.InsertKey(node);
    end
else
    % initial bounds
    I = eye(n);
    L = -inf;
    ub0 = inf;
    for i = 1:n
        p0 = replace(p, x, I(:, i));
        ub0 = min(ub0, p0);
    end
    U = ub0;
    
    % initial nodes
    H = MinHeap_BnB(opts.maxit*2);
    m = 2^(n-1);
    e = [2*(dec2bin(m-1:-1:0) - '0') - 1, ones(m,1)]; % half of binary vectors
    for i = 1:m
        V = I .* e(i, :);
        lb = get_lower_bound(x, p, n, d, V, opts);
        L = max(L, lb);
        center = mean(V, 2);
        ub = min(ub0, replace(p, x, center));
        U = min(U, ub);
        node = Node([], [], lb, ub, V);
        H.InsertKey(node);
    end
end

if opts.time
    time1 = toc(time0);
    out.tvec = [time1];
end

out.bnb = [];
out.lb_vec = [];
out.ub_vec = [U];
step = 1;

while (step < opts.maxit) && ~H.IsEmpty()
    node = H.ExtractMin();
    L = node.lb;
    out.lb_vec = [out.lb_vec L];
    if opts.verbose
        disp(['Step ', num2str(step), ': global lb = ', num2str(L)])
    end
    if U - L < opts.eps*(1 + abs(L) + abs(U))
        break
    end

    if opts.time
        time1 = toc(time0);
        out.tvec = [out.tvec time1];
    end

    [i, j, w] = partition(node);
    node_new1 = new_node(x, p, n, d, i, w, node, opts);
    node_new2 = new_node(x, p, n, d, j, w, node, opts);

    % upper bound update
    U = min([U, node_new1.ub, node_new2.ub]);
    out.ub_vec = [out.ub_vec U];
    if opts.verbose
        disp(['Step ', num2str(step), ': global ub = ', num2str(U)])
    end

    % pruning    
    if node_new1.lb < U
        H.InsertKey(node_new1);
        if opts.verbose
            disp(['Step ', num2str(step), ': lb = ', num2str(node_new1.lb), ', ub = ', num2str(node_new1.ub)])
        end
    else
        if opts.verbose
            disp(['Step ', num2str(step), ': lb = ', num2str(node_new1.lb), ' (pruned)'])
        end
    end
    if node_new2.lb < U
        H.InsertKey(node_new2);
        if opts.verbose
            disp(['Step ', num2str(step), ': lb = ', num2str(node_new2.lb), ', ub = ', num2str(node_new2.ub)])
        end
    else
        if opts.verbose
            disp(['Step ', num2str(step), ': lb = ', num2str(node_new2.lb), ' (pruned)'])
        end
    end
    
    step = step + 1;
end
t = L;
out.iter = step;
out.bnb = H;
end

function [i, j, w] = partition(node)
V = node.x;
cross_prod = V'*V;
min_prod = min(cross_prod, [], 'all');
[i, j] = find(cross_prod == min_prod, 1);
w = (V(:, i) + V(:, j))/2;
w = w / norm(w, 2);
end

function node_new = new_node(x, p, n, d, i, w, node, opts)
V = node.x;
V(:, i) = w;
lb = get_lower_bound(x, p, n, d, V, opts);
center = mean(V, 2);
ub = min(replace(p, x, w), replace(p, x, center));
node_new = Node([], [], lb, ub, V);
end

function lb = get_lower_bound(x, p, n, d, V, opts)
if opts.verbose == 2
    options = sdpsettings('solver', opts.solver, 'verbose', 1);
else
    options = sdpsettings('solver', opts.solver, 'verbose', 0);
end

gamma = sdpvar();
y = sdpvar(n, 1);
if opts.method == 0
    z = V * (y.^2);
    p_new = replace(p, x, z);
    reg_new = sum(z.^2)^(d/2);
    F = [sos(p_new - gamma*reg_new)];
    [sol, ~, ~] = solvesos(F, -gamma, options, [gamma]);
    lb = value(gamma);
elseif opts.method == 1
    z = V * y;
    p_new = replace(p, x, z);
    reg_new = sum(z.^2)^(d/2);
    F = [coefficients(p_new - gamma*reg_new, y) >= 0];
    sol = optimize(F, -gamma, options);
    if sol.problem == 1
        lb = -inf;
    else
        lb = value(gamma);
    end
end
end