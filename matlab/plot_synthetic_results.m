% analyse_synthetic_new3

DATASET = 1;

DECAYIND = 3;
DELAYIND = 5;

DELAYS = [0, 10, 20, 30];
HALFLIFES = [2 4 8 16 32 64];

boxmat = zeros(2000, 24, 2);

figure(1);
for k=1:length(DELAYS),
  for l=1:6,
    %subplot(6, 4, k+4*(l-1));
    %my_x = 2*(k+4*(l-1));
    my_x = 2*(l + 6*(k-1));
    v = sigmoidabTransform(squeeze(samples{l,k,DATASET}(:, DELAYIND, :)), 'atox', settings{DELAYIND});
    boxmat(:, my_x/2, 1) = v(:);
    [f, xi] = ksdensity(v(:));
    f_scaled = 0.95 * f / max(f);
    plot(my_x - f_scaled, xi)
    hold on
    plot(my_x + f_scaled, xi)
    %plot(xi, f);
    %axis([-5 80 0 1.1*max(f)])
    %hold on
    %plot(DELAYS(k)*[1, 1], [0, 1.1*max(f)], 'r')
  end
  plot(12*(k-1)+[1, 13], DELAYS(k)*[1, 1], 'k', 'LineWidth', 2)
end
axis([0 50 -5 80])
ylabel('\Delta (min)')
set(gca, 'XTick', 2:2:48)
set(gca, 'XTickLabel', [1:6, 1:6, 1:6, 1:6])
hold off
print -deps2 synthetic_delays_violins_2014-05-05


figure(2);
for l=1:6,
  for k=1:length(DELAYS),
    %subplot(6, 4, k+4*(l-1));
    my_x = 2*(k + 4*(l-1));
    v = 1./sigmoidabTransform(squeeze(samples{l,k,DATASET}(:, DECAYIND, :)), 'atox', settings{DECAYIND});
    boxmat(:, my_x/2, 2) = v(:);
    [f, xi] = ksdensity(v(:), 'bandwidth', 5);
    f_scaled = 0.95 * f / max(f);
    plot(my_x - f_scaled, xi)
    hold on
    plot(my_x + f_scaled, xi)
    %plot(xi, f);
    %axis([-5 80 0 1.1*max(f)])
    %hold on
    %plot(HALFLIFES(l)*[1, 1], [0, 1.1*max(f)], 'r')
    %hold off
  end
  plot(8*(l-1)+[1, 9], HALFLIFES(l)*[1, 1], 'k', 'LineWidth', 2)
end
axis([0 50 -5 80])
ylabel('t_{1/2} (min)')
set(gca, 'XTick', 2:2:48)
set(gca, 'XTickLabel', [1:4, 1:4, 1:4, 1:4, 1:4, 1:4])
hold off
print -deps2 synthetic_halflives_violins_2014-05-05

figure(3);
boxplot(boxmat(:, :, 1), 'symbol', '');
ylabel('\Delta (min)')
axis([0.5, 24.5, -2, 80])
hold on
for k=1:length(DELAYS),
  plot(6*(k-1)+[0.5, 6.5], DELAYS(k)*[1, 1], 'k', 'LineWidth', 2)
end
print -depsc2 synthetic_delays_2014-05-05

figure(4);
boxplot(boxmat(:, :, 2), 'symbol', '');
ylabel('t_{1/2} (min)')
axis([0.5, 24.5, -2, 80])
hold on
for k=1:length(HALFLIFES),
  plot(4*(k-1)+[0.5, 4.5], HALFLIFES(k)*[1, 1], 'k', 'LineWidth', 2)
end
print -depsc2 synthetic_halflives_2014-05-05
