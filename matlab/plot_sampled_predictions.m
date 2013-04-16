function plot_sampled_predictions(m, HMCsamples, gene_name, NODELAYHIST, SQRTTIME),

if nargin < 4,
  NODELAYHIST=0;
end

if nargin < 5,
  SQRTTIME=1;
end

shadecols = {[1.0 0.9 0.9], [0.9 1.0 0.9]};
linecols = {'r', 'g'};
FONTSIZE=9;
LINEWIDTH=2;

titles = {'Pol2 (input)', 'mRNA'};

delayI = 5;
settings = gpnddisimExtractParamTransformSettings(m);

if SQRTTIME,
  t_pred = (((0:100)/100*sqrt(1280)).^2 + 300)';
  t_tick = m.t{1} - min(m.t{1});
  t_plot = sqrt(t_pred - min(t_pred));
  t_tickplot = sqrt(t_tick);
else
  t_pred = [linspace(0, 157.5, 64), exp(linspace(log(160), log(1280), 61))]';
  t_tick = m.t{1} - min(m.t{1});
  t_plot = (1:length(t_pred))';
  t_tickplot = zeros(size(t_tick));
  for k=1:length(t_tick)
    t_tickplot(k) = find(abs(t_pred - t_tick(k)) < 1);
  end
  t_pred = t_pred + 300;
end

r = gpnddisimSamplePredictions(m, HMCsamples,t_pred, 500);
p = prctile(r, [2.5 97.5]);
mu = mean(r);

for k=1:2,
  if NODELAYHIST,
    subplot(2, 1, k);
  else
    subplot(2, 3, 3*k-2:3*k-1);
  end
  I = (1:length(t_pred)) + (k-1)*length(t_pred);
  J = (1:length(t_tick)) + (k-1)*length(t_tick);
  
  mycol = shadecols{k};
  h=fill([t_plot; t_plot(end:-1:1)], [p(1,I)'; flipud(p(2,I)')],mycol);
  set(h,'EdgeColor',mycol);
  %set(h,'FaceAlpha',0.5);
  %set(h,'EdgeAlpha',0.5);
  hold on
  plot(t_plot, mean(r(:, I)), linecols{k}, 'LineWidth', LINEWIDTH);
  switch k,
   case 1,
    plot(t_tickplot, m.y(J), 'bo')
   case 2,
    errorbar(t_tickplot, m.y(J), 2*sqrt(m.kern.comp{2}.comp{2}.fixedvariance), 'bo')
  end
  hold off
  V = axis;
  switch k,
    case 1,
      axis([0 max(t_plot) -1, ceil(max(m.y(J)))+1]);
    case 2,
      axis([0 max(t_plot) 0, ceil(max(m.y(J)))+1]);
  end
  %axis([0 max(t_plot) V(3:4)]);
  set(gca, 'XTick', t_tickplot);
  set(gca, 'XTickLabel', t_tick)
  set(gca, 'FontSize', FONTSIZE);
  ylabel(titles{k});
end

if ~NODELAYHIST,
  subplot(2, 3, [3 6]);
  vals = sigmoidabTransform(HMCsamples(:, delayI), 'atox', settings{delayI});
  [n, x] = hist(vals, 5:10:295);
  bar(x,n/sum(n),settings{delayI},'hist');
  if x(max(find(n))) < 100,
    V = axis;
    axis([0 100 V(3:4)]);
  end
  set(gca, 'FontSize', FONTSIZE);
  xlabel('Delay (min)')
end

if nargin > 2,
  if NODELAYHIST,
    subplot(2, 1, 1);
    title(gene_name);
    subplot(2, 1, 2);
  else
    subplot(2, 3, 1:2);
    title(gene_name);
    subplot(2, 3, [3 6]);
  end
end
