function plot_data(pol, rna, gene),

times = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280];

subplot(2, 1, 1);
plot(sqrt(times), pol, 'x-')
title(gene)
set(gca, 'XTick', sqrt(times))
set(gca, 'XTickLabel', times)
subplot(2, 1, 2);
plot(sqrt(times), rna, 'x-')
set(gca, 'XTick', sqrt(times))
set(gca, 'XTickLabel', times)
