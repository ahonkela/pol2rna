function plot_sampled_predictions(m, HMCsamples, gene_name, NODELAYHIST, SQRTTIME, premrnadata, SHORTHIST),

if nargin < 4,
  NODELAYHIST=0;
end

if nargin < 5,
  SQRTTIME=1;
end

if nargin < 6,
  PLOTPREMRNA=0;
else
  if isempty(premrnadata),
    PLOTPREMRNA=0;
  else
    PLOTPREMRNA=1;
  end
end

if nargin < 7,
  SHORTHIST = 0;
end

if PLOTPREMRNA,
  PLOTROWS = 3;
else
  PLOTROWS = 2;
end

if NODELAYHIST,
  PLOTCOLS = 1;
else
  PLOTCOLS = 3;
end

shadecols = {[1.0 0.9 0.9], [0.9 1.0 0.9]};
linecols = {'r', 'g', 'b'};
FONTSIZE=6;
LINEWIDTH=1.5;

titles = {'Pol II', 'mRNA (FPKM)', 'pre-mRNA'};

delayI = 5;
settings = gpnddisimExtractParamTransformSettings(m);

t_len = m.t{1}(end) - m.t{1}(1);

switch SQRTTIME,
  case 2,
    t_pred = linspace(m.t{1}(1), m.t{1}(end), 101)';
    t_tick = m.t{1} - min(m.t{1});
    t_plot = t_pred - min(t_pred);
    t_tickplot = t_tick;
    t_ticklabels = t_tick;
  case 1,
    t_pred = (((0:100)/100*sqrt(t_len)).^2 + 300)';
    t_tick = m.t{1} - min(m.t{1});
    t_plot = sqrt(t_pred - min(t_pred));
    t_tickplot = sqrt(t_tick);
    t_ticklabels = t_tick;
  case 0,
    t_pred = [linspace(0, 157.5, 64), exp(linspace(log(160), log(1280), 61))]';
    t_tick = m.t{1} - min(m.t{1});
    t_plot = (1:length(t_pred))';
    t_tickplot = zeros(size(t_tick));
    t_ticklabels = cell(size(t_tick));
    for k=1:length(t_tick)
      t_tickplot(k) = find(abs(t_pred - t_tick(k)) < 1);
      t_ticklabels{k} = sprintf('%d', t_tick(k));
    end
    t_ticklabels{2} = '';
    t_ticklabels{3} = '';
    t_ticklabels{4} = '';
    t_pred = t_pred + 300;
end

r = gpnddisimSamplePredictions(m, HMCsamples,t_pred, 500);
p = prctile(r, [2.5 97.5]);
mu = mean(r);

for k=1:2,
  if NODELAYHIST,
    if PLOTPREMRNA,
      subplot(PLOTROWS, PLOTCOLS, 2*k-1);
    else
      subplot(PLOTROWS, PLOTCOLS, k);
    end
  else
    if PLOTPREMRNA,
      subplot(PLOTROWS, PLOTCOLS, 6*k-5:6*k-4);
    else
      subplot(PLOTROWS, PLOTCOLS, 3*k-2:3*k-1);
    end
  end
  I = (1:length(t_pred)) + (k-1)*length(t_pred);
  J = (1:length(t_tick)) + (k-1)*length(t_tick);

  % plot bounds
  switch k,
    case 1,
      minbound = 0;
      maxbound = 1.2*(max([max(m.y(J)), max(mu(I))]) + 1e-3);
    case 2,
      minbound = 0;
      maxbound = ceil(1.1 * max(m.y(J) ...
                                + 2*sqrt(m.kern.comp{2}.comp{2}.fixedvariance)));
  end

  % clip error bars to plot bounds
  p(1,I(p(1,I) < minbound)) = minbound;
  p(2,I(p(2,I) > maxbound)) = maxbound;
  
  mycol = shadecols{k};
  h=fill([t_plot; t_plot(end:-1:1)], [p(1,I)'; flipud(p(2,I)')],mycol);
  set(h,'EdgeColor',mycol);
  %set(h,'FaceAlpha',0.5);
  %set(h,'EdgeAlpha',0.5);
  hold on
  plot(t_plot, mean(r(:, I)), linecols{k}, 'LineWidth', LINEWIDTH);
  switch k,
   case 1,
    plot(t_tickplot, m.y(J), 'bx')
   case 2,
    errorbar(t_tickplot, m.y(J), 2*sqrt(m.kern.comp{2}.comp{2}.fixedvariance), 'bx')
  end
  hold off
  axis([0, max(t_plot), minbound, maxbound]);
  set(gca, 'XTick', t_tickplot);
  if k==2,
    set(gca, 'XTickLabel', t_ticklabels)
  else
    set(gca, 'XTickLabel', []);
  end
  set(gca, 'FontSize', FONTSIZE);
  v = axis;
  text(v(1)-0.22*(v(2)-v(1)), mean(v(3:4)), titles{k}, 'HorizontalAlignment','center', 'Rotation', 90, 'FontSize', FONTSIZE)
  if k==2,
    v = axis;
    text(mean(v(1:2)), v(3)-0.26*(v(4)-v(3)), 't (min)', 'HorizontalAlignment','center', 'FontSize', FONTSIZE);
  end
end

if PLOTPREMRNA,
  if NODELAYHIST,
    subplot(PLOTROWS, PLOTCOLS, 2);
  else
    subplot(PLOTROWS, PLOTCOLS, 4:5);
  end
  plot(t_tickplot, premrnadata, [linecols{3} 'x-'])
  axis([0, max(t_plot), 0, 1.2*(max(premrnadata)+1e-3)])
  set(gca, 'XTick', t_tickplot);
  set(gca, 'XTickLabel', [])
  set(gca, 'FontSize', FONTSIZE);
  v = axis;
  text(v(1)-0.22*(v(2)-v(1)), mean(v(3:4)), titles{3}, 'HorizontalAlignment','center', 'Rotation', 90, 'FontSize', FONTSIZE)
end

if ~NODELAYHIST,
  subplot(PLOTROWS, PLOTCOLS, PLOTCOLS*(1:PLOTROWS));
  vals = sigmoidabTransform(HMCsamples(:, delayI), 'atox', settings{delayI});
  [n, x] = hist(vals, 5:10:295);
  if x(max(find(n))) < 100 | SHORTHIST,
    %V = axis;
    %axis([0 100 V(3:4)]);
    bar(x(1:10),n(1:10)/sum(n),[0, 100],'hist');
    axis([0 100 0 1]);
    xmax = 100;
  else
    bar(x,n/sum(n),settings{delayI},'hist');
    axis([0 300 0 1]);
    xmax = 300;
  end
  set(gca, 'FontSize', FONTSIZE);
  if xmax == 100,
    set(gca, 'XTick', 0:50:xmax);
    set(gca, 'XTickLabel', 0:50:xmax);
  else
    set(gca, 'XTick', 0:100:xmax);
    set(gca, 'XTickLabel', 0:100:xmax);
  end
  v = axis;
  text(mean(v(1:2)), v(3)-0.10*(v(4)-v(3)), 'Delay (min)', 'HorizontalAlignment','center', 'FontSize', FONTSIZE);
end

if nargin > 2,
  if NODELAYHIST,
    subplot(PLOTROWS, PLOTCOLS, 1);
    title(gene_name);
    subplot(PLOTROWS, PLOTCOLS, PLOTROWS);
  else
    subplot(PLOTROWS, PLOTCOLS, 1:2);
    title(gene_name);
    subplot(PLOTROWS, PLOTCOLS, PLOTCOLS*(1:PLOTROWS));
  end
end
