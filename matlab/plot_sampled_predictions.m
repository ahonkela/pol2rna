function plot_sampled_predictions(m, HMCsamples, gene_name),

shadecols = {[1.0 0.7 0.7], [0.7 1.0 0.7]};
linecols = {'r', 'g'};
FONTSIZE=9;
LINEWIDTH=2;

titles = {'Pol2 (input)', 'mRNA'};

delayI = 5;
settings = gpnddisimExtractParamTransformSettings(m);

t_pred = (((0:100)/100*sqrt(1280)).^2 + 300)';
t_tick = m.t{1} - min(m.t{1});
t_plot = sqrt(t_pred - min(t_pred));
t_tickplot = sqrt(t_tick);

r = gpnddisimSamplePredictions(m, HMCsamples,t_pred, 500);
p = prctile(r, [2.5 97.5]);
mu = mean(r);

for k=1:2,
  subplot(3, 2, 2*k-1:2*k);
  I = (1:length(t_pred)) + (k-1)*length(t_pred);
  J = (1:length(t_tick)) + (k-1)*length(t_tick);
  
  mycol = shadecols{k};
  h=fill([t_plot; t_plot(end:-1:1)], [p(1,I)'; flipud(p(2,I)')],mycol);
  set(h,'EdgeColor',mycol);
  set(h,'FaceAlpha',0.5);
  set(h,'EdgeAlpha',0.5);
  hold on
  plot(sqrt(t_pred - min(t_pred)), mean(r(:, I)), linecols{k}, 'LineWidth', LINEWIDTH);
  switch k,
   case 1,
    plot(t_tickplot, m.y(J), 'bo-')
   case 2,
    errorbar(t_tickplot, m.y(J), 2*sqrt(m.kern.comp{2}.comp{2}.fixedvariance), 'bo-')
  end
  hold off
  V = axis;
  axis([0 max(t_plot) V(3:4)]);
  set(gca, 'XTick', t_tickplot);
  set(gca, 'XTickLabel', t_tick)
  set(gca, 'FontSize', FONTSIZE);
  ylabel(titles{k});
end

subplot(3, 2, 5);
vals = sigmoidabTransform(HMCsamples(:, delayI), 'atox', settings{delayI});
[n, x] = hist(vals, 5:10:295);
bar(x,n/sum(n),settings{delayI},'hist');
set(gca, 'FontSize', FONTSIZE);
xlabel('Delay (min)')

if nargin > 2,
  subplot(3, 2, 1:2);
  title(gene_name);
  subplot(3, 2, 5);
end
