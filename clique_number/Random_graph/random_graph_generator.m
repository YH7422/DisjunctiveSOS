n = 75;
p = 0.5;
m = 5;
output_dir = '.';

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for k = 1:m
    upper = triu(rand(n) < p, 1);
    A = upper + upper';
    fname = "Erdos_n" + num2str(n) + "_p" + num2str(p) + "_test" + num2str(k) + ".mat";
    save(fullfile(output_dir, fname), 'A');
end