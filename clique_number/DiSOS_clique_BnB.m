% A branch-and-bound method aimimg to compute the upper bound of 
% the clique number of graph G by solving
% min t s.t. t*(I+\bar{A}) - J \in C_n

function [t, out] = DiSOS_clique_BnB(A, opts)
if ~isfield(opts, 'max_node'); opts.max_node = 1000; end
if ~isfield(opts, 'maxit'); opts.maxit = 5; end
if ~isfield(opts, 'eps'); opts.eps = 1e-5; end
if ~isfield(opts, 'solver'); opts.solver = 'MOSEK'; end
if ~isfield(opts, 'dd'); opts.dd = 0; end
% dd = 1 means using diagonally dominant matrices to 
% approximate the PSD cone
if ~isfield(opts, 'verbose'); opts.verbose = 1; end
% verbose = 0 no output
% verbose = 1 basic output
% verbose = 2 detailed output
if ~isfield(opts, 'time'); opts.time = 0; end
if ~isfield(opts, 'seed'); opts.seed = 0; end
if ~isfield(opts, 'name'); opts.name = ""; end

% fix the seed to reproduce the result
ss = RandStream('mt19937ar', 'Seed', opts.seed); 
RandStream. setGlobalStream(ss);

% CVX settings
cvx_solver MOSEK
if opts.verbose <= 1
    cvx_quiet true
else
    cvx_quiet false
end

n = size(A, 1);
I = eye(n);
out = struct();
out.bnb = [];
out.lb_vec = [];
out.ub_vec = [];
if opts.time
    time0 = tic;
    out.tvec = [];
end

% initial node
[t, i0, j0, delta] = get_lower_bound(I, A, opts);
f0 = max(diag(A));
lb = -t;
ub = -1 / (1-f0);
node = Node([], [], lb, ub, I);
node.cutting_edge = [i0, j0];

H = MinHeap_BnB(opts.max_node);
H.InsertKey(node);
L = lb;
U = ub;
step = 1;

while (step < opts.max_node) && ~H.IsEmpty()
    node = H.ExtractMin();
    L = node.lb;
    out.lb_vec = [out.lb_vec L];
    if opts.verbose
        disp(['Step ', num2str(step), ': global ub = ', num2str(-L)])
    end
    if floor(-L) == ceil(-U)
        break
    end

    if opts.time
        time1 = toc(time0);
        out.tvec = [out.tvec time1];
    end

    i = node.cutting_edge(1);
    if i == -1
       % The lower bound is exact
        if opts.verbose
            disp(['Step ', num2str(step), ': the lower bound is exact'])
        end
        continue
    end
    j = node.cutting_edge(2);
    x = node.x;
    w = (x(i, :) + x(j, :))/2;
    node_new1 = new_node(A, i, w, node, opts);
    node_new2 = new_node(A, j, w, node, opts);
    
    % upper bound update
    U = min([U, node_new1.ub, node_new2.ub]);
    out.ub_vec = [out.ub_vec U];
    if opts.verbose
        disp(['Step ', num2str(step), ': global lb = ', num2str(-U)])
    end

    % pruning    
    if node_new1.lb < U
        H.InsertKey(node_new1);
        if opts.verbose
            disp(['Step ', num2str(step), ': lb = ', num2str(-node_new1.ub), ', ub = ', num2str(-node_new1.lb)])
        end
    else
        if opts.verbose
            disp(['Step ', num2str(step), ': pruned'])
        end
    end
    if node_new2.lb < U
        H.InsertKey(node_new2);
        if opts.verbose
            disp(['Step ', num2str(step), ': lb = ', num2str(-node_new2.ub), ', ub = ', num2str(-node_new2.lb)])
        end
    else
        if opts.verbose
            disp(['Step ', num2str(step), ': pruned'])
        end
    end

    step = step + 1;
end

if H.IsEmpty()
    out.lb_vec = [out.lb_vec U];
end
t = L;
out.iter = step;
out.bnb = H;
end

function node_new = new_node(A, i, w, node, opts)
n = size(A, 1);

% get lower bound
x = node.x;
x(i, :) = w;
[t, i0, j0, delta] = get_lower_bound(x, A, opts);
lb = -t;

% get upper bound
if i0 == -1
    % The lower bound is exact
    ub = lb;
else
    alpha = zeros(n, 1);
    alpha(i) = 1;
    fx = w*A*w';
    % use projected gradient descent to maximize a'*G*a on the unit simplex
    G = x*A*x';
    ss = delta * 1e-1;
    for i = 1 : opts.maxit
        alpha = proj_prob_vec(alpha + ss*G*alpha, 0, 1);
        fx = max(fx, alpha'*G*alpha);
    end
    ub = -1 / (1-fx);
end

node_new = Node([], [], lb, ub, x);
node_new.cutting_edge = [i0, j0];
end

function [t, i, j, delta] = get_lower_bound(x, A, opts)
n = size(A, 1);
I = eye(n);
J = ones(n);
Ac = J - I - A;
if opts.dd
    cvx_begin
        variables t N(n, n)
        minimize(t)
        D = x*(t*(I+Ac) - J)*x' - N;
        2 * diag(D) >= sum(abs(D), 2);
        N >= 0;
    cvx_end
else
    cvx_begin
        variables t N(n, n)
        minimize(t)
        x*(t*(I+Ac) - J)*x' - N == semidefinite(n);
        N >= 0;
    cvx_end
end

[i_vec, j_vec] = find(abs(N) <= opts.eps*max(max(N, [], 'all'), 1));
dist = pdist(x);
delta = max(dist);
D = squareform(dist);
max_dist_active = -inf;
i = -1;
j = -1;
for k = 1:length(i_vec)
    i0 = i_vec(k);
    j0 = j_vec(k);
    if i0 >= j0
        continue
    end
    if max_dist_active < D(i0,j0)
        max_dist_active = D(i0,j0);
        i = i0;
        j = j0;
    end
end
end

function x = proj_prob_vec(y, l, r)
n = size(y, 1);
x = zeros(n, 1);
tot = sum(l);
lambda = [l-y; r-y];
[~, idx] = sort(lambda);
lambda = lambda(idx);
active = 1;
for i = 2 : n*2
    tot = tot + active * (lambda(i) - lambda(i-1));
    if tot >= 1
        lam = (1-tot)/active + lambda(i);
        x = min(max(y + lam, l), r);
        return
    elseif idx(i) <= n
        active = active + 1;
    else
        active = active - 1;
    end
end
if all(x == 0)
    error('Incorrect projection.')
end
end