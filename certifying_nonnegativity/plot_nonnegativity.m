name = 'Robinson2';
init = 1;
n = 4;
filename = [char(name), '_bnb_sd_init', num2str(init)];
load(filename);
if init == 0
    m = 2^(n-1)-1;
else
    m = n;
end

figure('Visible', 'off')
p1 = plot(m+(1:length(out.lb_vec)), out.lb_vec, 'DisplayName', 'lower bound', 'LineWidth', 1);
hold on
p2 = plot(m+(1:length(out.ub_vec)), out.ub_vec, '-.', 'DisplayName', 'upper bound', 'LineWidth', 1);
% yline(0)
xlabel('Number of regions')
ylabel('Objective value')
legend([p1, p2], 'Location', 'southeast')
ylim([-inf, max([max(out.ub_vec)*1.5, abs(out.lb_vec(1))*0.1, 1e-4])])
ax = gca;
ax.XTick = unique(round(ax.XTick));
ax.FontSize = 12;
exportgraphics(gcf,[filename, '.pdf'], 'ContentType','vector')
