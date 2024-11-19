n = 3;
x = sdpvar(n, 1);
opts = struct();
opts.init = 1;
opts.verbose = 0;
opts.maxit = 300;
opts.dq = 2;
for m = 1:6
    d = 2*(2*m+1);
    p = x(1)^(2*m+1) * x(3)^(2*m+1) + (x(2)^2 * x(3)^(2*m-1) - x(1)^(2*m+1) - x(1)*x(3)^(2*m))^2;
    [t, out] = DiSOS_BnB_SD(x, p, n, d, opts);
end