test = 4;
filename = ['copositive', num2str(test), '_bnb'];
load(filename);

figure('Visible', 'off')
p1 = plot((2:length(out.lb_vec)), out.lb_vec(2:end), 'Displayname', 'lower bound', 'LineWidth', 1);
hold on
p2 = plot((2:length(out.ub_vec)), out.ub_vec(2:end), '-.', 'Displayname', 'upper bound', 'LineWidth', 1);
xlabel('Number of regions');
ylabel('Objective value')
legend([p1, p2], 'Location', 'northeast')
ax = gca;
ax.FontSize = 12;
exportgraphics(gcf,[filename, '.pdf'], 'ContentType','vector')
