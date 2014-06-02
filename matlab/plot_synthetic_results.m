% analyse_synthetic_new3

FONTSIZE = 6;
LABELSHIFT = 10;

DATASET = 1;

DECAYIND = 3;
DELAYIND = 5;

DELAYS = [0, 10, 20, 30];
HALFLIFES = [2 4 8 16 32 64];

boxmat = zeros(2000, 24, 2);

for k=1:length(DELAYS),
  for l=1:6,
    my_x = 2*(l + 6*(k-1));
    v = sigmoidabTransform(squeeze(samples{l,k,DATASET}(:, DELAYIND, :)), 'atox', settings{DELAYIND});
    boxmat(:, my_x/2, 1) = v(:);
  end
end

for l=1:6,
  for k=1:length(DELAYS),
    %subplot(6, 4, k+4*(l-1));
    my_x = 2*(k + 4*(l-1));
    v = 1./sigmoidabTransform(squeeze(samples{l,k,DATASET}(:, DECAYIND, :)), 'atox', settings{DECAYIND});
    boxmat(:, my_x/2, 2) = v(:);
  end
end

figure(1);
set(gca, 'FontSize', FONTSIZE)
%boxplot(boxmat(:, :, 1), 'symbol', '', 'labels', repmat([2 4 8 16 32 64], [1 4]));
bplot(boxmat(:, :, 1), 'nomean', 'linewidth', 1);
set(gca, 'FontSize', FONTSIZE)
ylabel('\Delta (min)')
set(gca, 'XTick', 1:24);
set(gca, 'XTickLabel', repmat([2 4 8 16 32 64], [1 4]));
axis([0.5, 24.5, -1, 80])
hold on
for k=1:length(DELAYS),
  plot(6*(k-1)+[0.5, 6.5], DELAYS(k)*[1, 1], 'k', 'LineWidth', 2)
end
hold off
set(gca, 'FontSize', FONTSIZE)
%set(findobj(gca,'Type','text'),'FontSize',FONTSIZE);
%set(findobj(h,'Type','text'),'VerticalAlignment','top');
%set(gca, 'XTickLabel', repmat([2 4 8 16 32 64], [1 4]))
xlabel('t_{1/2} (min)')
set(gcf, 'PaperUnits', 'centimeters')
set(gcf,'PaperPositionMode','auto')
set(gcf, 'PaperPosition', [0 0 8.7, 7.0])
print -depsc2 synthetic_delays

figure(2);
set(gca, 'FontSize', FONTSIZE)
bplot(boxmat(:, :, 2), 'nomean', 'linewidth', 1);
set(gca, 'FontSize', FONTSIZE)
ylabel('t_{1/2} (min)')
set(gca, 'XTick', 1:24);
set(gca, 'XTickLabel', repmat([0 10 20 30], [1 6]));
axis([0.5, 24.5, -1, 80])
hold on
for k=1:length(HALFLIFES),
  plot(4*(k-1)+[0.5, 4.5], HALFLIFES(k)*[1, 1], 'k', 'LineWidth', 2)
end
hold off
set(gca, 'FontSize', FONTSIZE)
%set(findobj(gca,'Type','text'),'FontSize',FONTSIZE);
%set(gca, 'XTickLabel', repmat([0 10 20 30], [1 6]))
xlabel('\Delta (min)')
set(gcf, 'PaperUnits', 'centimeters')
set(gcf,'PaperPositionMode','auto')
set(gcf, 'PaperPosition', [0 0 8.7, 7.0])
print -depsc2 synthetic_halflives
