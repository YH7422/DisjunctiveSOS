% A branch-and-bound method aiming to compute the upper bound of 
% the clique number of graph G by solving
% min t s.t. t*(I+\bar{A}) - J \in C_n

function [t, out] = DiSOS_clique_BnB(A, opts)
if ~isfield(opts, 'max_node'); opts.max_node = 1000; end
if ~isfield(opts, 'maxit'); opts.maxit = 5; end
if ~isfield(opts, 'alpha'); opts.alpha = 0.5; end
if ~isfield(opts, 'eps'); opts.eps = 1e-5; end
if ~isfield(opts, 'dd'); opts.dd = 0; end
if ~isfield(opts, 'verbose'); opts.verbose = 1; end
if ~isfield(opts, 'time'); opts.time = 0; end
if ~isfield(opts, 'seed'); opts.seed = 0; end
if ~isfield(opts, 'name'); opts.name = ""; end

ss = RandStream('mt19937ar', 'Seed', opts.seed); 
RandStream.setGlobalStream(ss);

n = size(A, 1);
Ie = eye(n);
Ac = ones(n) - Ie - A;
B = Ie + Ac;
out = struct();
out.bnb = [];
out.lb_vec = [];
out.ub_vec = [];
if opts.time
    time0 = tic;
    out.tvec = [];
end

% Build a reusable MOSEK problem template.
% The SDP solved at each node is:
%   min  t
%   s.t. S + N = t*M1 - J,   S psd,  N >= 0 elementwise,
% where M1 = x*B*x' changes per node and J = x*ones(n)*x' = ones(n) is fixed.
%
% Scalar variables: x(1) = t (free), x(2:m+1) = upper triangle of N (nonneg).
% Bar variable:     barX (n x n PSD matrix S).
% Constraints (m equalities, one per upper-triangle entry (i,j), i<=j):
%   barX_{ij} + N_{ij} - t*M1_{ij} = -1.
%
% Per solve we only update the first column of prob.a.
msk = setup_mosek_prob(n);

% initial node
[t_val, N_val, i0, j0, delta] = get_lower_bound(Ie, B, msk, opts);
f0 = max(diag(A));
lb = -t_val;
ub = -1 / (1 - f0);
node = Node([], [], lb, ub, Ie);
node.cutting_edge = [i0, j0];

H = MinHeap_BnB(opts.max_node);
H.InsertKey(node);
L = lb;
U = ub;
out.ub_vec = [U];
step = 1;

while (step < opts.max_node) && ~H.IsEmpty()
    node = H.ExtractMin();
    L = node.lb;
    out.lb_vec = [out.lb_vec L];
    if opts.verbose
        disp(['Step ', num2str(step), ': global ub = ', num2str(-L)])
    end
    if floor(-L) == ceil(-U - 1e-4)
        break
    end

    if opts.time
        time1 = toc(time0);
        out.tvec = [out.tvec time1];
    end

    i = node.cutting_edge(1);
    if i == -1
        if opts.verbose
            disp(['Step ', num2str(step), ': the lower bound is exact'])
        end
        continue
    end
    j = node.cutting_edge(2);
    x = node.x;
    w = (x(i, :) + x(j, :)) / 2;
    node_new1 = new_node(A, B, i, w, node, opts, msk);
    node_new2 = new_node(A, B, j, w, node, opts, msk);
    
    U = min([U, node_new1.ub, node_new2.ub]);
    out.ub_vec = [out.ub_vec U];
    if opts.verbose
        disp(['Step ', num2str(step), ': global lb = ', num2str(-U)])
    end

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

%% MOSEK direct interface

function msk = setup_mosek_prob(n)
m = n * (n + 1) / 2;
numvar = 1 + m;

% Upper-triangle indices (i <= j), sorted by column then row.
[I_idx, J_idx] = find(triu(ones(n)));
[~, perm] = sortrows([J_idx, I_idx]);
I_idx = I_idx(perm);
J_idx = J_idx(perm);
linidx = sub2ind([n, n], I_idx, J_idx);

% Variable bounds: t free, N entries nonneg.
prob.blx = [-inf; zeros(m, 1)];
prob.bux = inf(numvar, 1);

% Objective: minimize t.
prob.c = [1; zeros(m, 1)];

% Fixed part of constraint matrix:
% each constraint k has coefficient 1 on variable x(1+k) = N_{I(k),J(k)}.
A_fixed = sparse((1:m)', (2:m+1)', ones(m, 1), m, numvar);

% RHS is fixed: M2 = x*J*x' = J always since rows of x sum to 1.
rhs = -ones(m, 1);
prob.blc = rhs;
prob.buc = rhs;

% Bar variable: one n x n PSD matrix S.
prob.bardim = n;

% Bar constraint coefficients (fixed across solves).
% To extract S_{ij} from barX in constraint k:
%   diagonal (i==j):     store (i, i) with value 1
%   off-diagonal (i<j):  store (j, i) with value 0.5
%     because MOSEK uses subk >= subl and counts off-diag twice.
offdiag = (I_idx ~= J_idx);
bara_val = ones(m, 1);
bara_val(offdiag) = 0.5;

prob.bara.subi = (1:m)';
prob.bara.subj = ones(m, 1);
prob.bara.subk = J_idx;
prob.bara.subl = I_idx;
prob.bara.val  = bara_val;

% Store metadata.
msk.prob    = prob;
msk.A_fixed = A_fixed;
msk.linidx  = linidx;
msk.m       = m;
msk.n       = n;
end

function [t_val, N_val] = solve_mosek(msk, M1)
n = msk.n;
m = msk.m;
linidx = msk.linidx;
prob = msk.prob;

% Update first column of A: coefficient of t is -M1_{ij}.
t_coeffs = -M1(linidx);
prob.a = msk.A_fixed + sparse((1:m)', ones(m, 1), t_coeffs, m, 1 + m);

% Solve.
[~, res] = mosekopt('minimize echo(0)', prob);

t_val = res.sol.itr.xx(1);

% Reconstruct N.
N_val = zeros(n);
N_val(linidx) = res.sol.itr.xx(2:end);
N_val = N_val + triu(N_val, 1)';
end

%% BnB helpers

function node_new = new_node(A, B, i, w, node, opts, msk)
n = size(A, 1);

x = node.x;
x(i, :) = w;
[t_val, N_val, i0, j0, delta] = get_lower_bound(x, B, msk, opts);
lb = -t_val;

if i0 == -1
    ub = lb;
else
    alpha = zeros(n, 1);
    alpha(i) = 1;
    fx = w * A * w';
    G = x * A * x';
%     alpha = ones(n, 1) / n;
%     fx = max(fx, alpha'*G*alpha);
    ss = delta * opts.alpha;
    for k = 1:opts.maxit
        alpha = proj_prob_vec(alpha + ss * G * alpha, 0, 1);
        fx = max(fx, alpha' * G * alpha);
    end
    ub = -1 / (1 - fx);
end

node_new = Node([], [], lb, ub, x);
node_new.cutting_edge = [i0, j0];
end

function [t_val, N_val, i, j, delta] = get_lower_bound(x, B, msk, opts)
M1 = x * B * x';

[t_val, N_val] = solve_mosek(msk, M1);

[i_vec, j_vec] = find(abs(N_val) <= opts.eps * max(max(N_val, [], 'all'), 1));
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
    if max_dist_active < D(i0, j0)
        max_dist_active = D(i0, j0);
        i = i0;
        j = j0;
    end
end
end

function x = proj_prob_vec(y, l, r)
n = size(y, 1);
x = zeros(n, 1);
tot = sum(l);
lambda = [l - y; r - y];
[~, idx] = sort(lambda);
lambda = lambda(idx);
active = 1;
for i = 2:n * 2
    tot = tot + active * (lambda(i) - lambda(i - 1));
    if tot >= 1
        lam = (1 - tot) / active + lambda(i);
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