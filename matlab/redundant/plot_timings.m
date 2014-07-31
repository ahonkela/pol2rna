obst = [0 5 10 20 40 80 160 320 640 1280];

[~, polmax] = max(mu(:, 1:101), [], 2);
[~, rnamax] = max(mu(:, 102:end), [], 2);
plot(polmax + 0.2*randn(size(polmax)) - 1, rnamax + 0.2*randn(size(rnamax)) - 1, '.')
axis([-1 101 -1 101])
xlabel('Pol2 peak time');
ylabel('RNA peak time');
set(gca, 'XTick', sqrt(obst) * 100 / sqrt(1280));
set(gca, 'XTickLabel', obst)
set(gca, 'YTick', sqrt(obst) * 100 / sqrt(1280));
set(gca, 'YTickLabel', obst)
set(gcf, 'PaperPosition', [0 0 12 12]);
