% analyse_synthetic_new3

DATASET = 1;

DECAYIND = 3;
DELAYIND = 5;

DELAYS = [0, 10, 20, 30];
HALFLIFES = [2 4 8 16 32 64];

figure(1);
for k=1:length(DELAYS),
  for l=1:6,
    subplot(6, 4, k+4*(l-1));
    v = sigmoidabTransform(squeeze(samples{l,k,DATASET}(:, DELAYIND, :)), 'atox', settings{DELAYIND});
    [f, xi] = ksdensity(v(:));
    plot(xi, f);
    axis([-5 80 0 1.1*max(f)])
    hold on
    plot(DELAYS(k)*[1, 1], [0, 1.1*max(f)], 'r')
    hold off
  end
end
print -depsc2 synthetic_delays_2014-05-05


figure(2);
for k=1:length(DELAYS),
  for l=1:6,
    subplot(6, 4, k+4*(l-1));
    v = 1./sigmoidabTransform(squeeze(samples{l,k,DATASET}(:, DECAYIND, :)), 'atox', settings{DECAYIND});
    [f, xi] = ksdensity(v(:), 'bandwidth', 5);
    plot(xi, f);
    axis([-5 80 0 1.1*max(f)])
    hold on
    plot(HALFLIFES(l)*[1, 1], [0, 1.1*max(f)], 'r')
    hold off
  end
end
print -depsc2 synthetic_halflives_2014-05-05
