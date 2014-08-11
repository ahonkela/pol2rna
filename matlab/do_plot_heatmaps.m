function I_pcomb = do_plot_heatmaps(dataVals1, dataVals2, name, Iextra, highlights),

if nargin < 5,
  highlights = [];
end

PRINTSTYLE='-painters';

sz = size(dataVals1);
dataVals1 = dataVals1 ./ repmat(max(dataVals1, [], 2), [1, sz(2)]);
dataVals2 = dataVals2 ./ repmat(max(dataVals2, [], 2), [1, sz(2)]);

% val1mean = mean(dataVals1 .* repmat(1:sz(2), [sz(1), 1]), 2);
% val2mean = mean(dataVals2 .* repmat(1:sz(2), [sz(1), 1]), 2);
% [~, I] = sort(val1mean + val2mean);
% plot_with_order(dataVals1, dataVals2, I);
% print('-depsc2', PRINTSTYLE, [name, '_heatmap_meansort.eps']);

% [~, I] = sort(val1mean);
% plot_with_order(dataVals1, dataVals2, I);
% print('-depsc2', PRINTSTYLE, [name, '_heatmap_pol2sort.eps']);

% [~, I] = sort(val2mean);
% plot_with_order(dataVals1, dataVals2, I);
% print('-depsc2', PRINTSTYLE, [name, '_heatmap_mrnasort.eps']);

[~, peak1] = max(dataVals1, [], 2);
[~, peak2] = max(dataVals2, [], 2);

% [~, I] = sort(peak1 + peak2);
% plot_with_order(dataVals1, dataVals2, I);
% print('-depsc2', PRINTSTYLE, [name, '_heatmap_pmeansort.eps']);

[~, I] = sort(1000*peak1 + peak2);
plot_with_order(dataVals1, dataVals2, I, highlights);
print('-depsc2', PRINTSTYLE, [name, '_heatmap_pcombsort.eps']);

I_pcomb = I;

% [~, I] = sort(peak1);
% plot_with_order(dataVals1, dataVals2, I);
% print('-depsc2', PRINTSTYLE, [name, '_heatmap_ppol2sort.eps']);

% [~, I] = sort(peak2);
% plot_with_order(dataVals1, dataVals2, I);
% print('-depsc2', PRINTSTYLE, [name, '_heatmap_pmrnasort.eps']);

% Y = pdist([dataVals1, dataVals2]);
% Z = linkage(Y, 'average');
% order = optimalleaforder(Z, Y);

% plot_with_order(dataVals1, dataVals2, order);
% print('-depsc2', PRINTSTYLE, [name, '_heatmap_hclustsort.eps']);

if nargin > 3 && ~isempty(Iextra),
  plot_with_order(dataVals1, dataVals2, Iextra, highlights);
  print('-depsc2', PRINTSTYLE, [name, '_heatmap_extrasort.eps']);
end



function plot_with_order(dataVals1, dataVals2, I, highlights)

FONTSIZE = 6;
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperPosition', [0, 0, 60/25.4, 100/25.4])

assert(size(dataVals1, 2) == 101);
t_tick = [0, 20, 80, 320, 640, 1280];
x_tick = sqrt(t_tick) * 100 / sqrt(1280) + 1;

colormap('hot')
subplot(1, 2, 1);
imagesc(dataVals1(I, :))
set(gca, 'FontSize', FONTSIZE);
if nargin > 3 && ~isempty(highlights),
  [J, JJ] = sort(I(highlights.I));
  labels = highlights.txt(JJ);
  set(gca, 'YTick', J);
  set(gca, 'YTickLabel', labels);
else
  set(gca, 'YTick', [])
end
set(gca, 'XTick', x_tick)
set(gca, 'XTickLabel', t_tick)
xlabel('t (min)')
title('Pol-II')
subplot(1, 2, 2);
imagesc(dataVals2(I, :))
set(gca, 'FontSize', FONTSIZE);
set(gca, 'YTick', [])
set(gca, 'XTick', x_tick)
set(gca, 'XTickLabel', t_tick)
xlabel('t (min)')
title('mRNA')
