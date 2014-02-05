function plot_sampled_predictions_premrna(m, HMCsamples, gene_name, NODELAYHIST, SQRTTIME, pol2data),

if nargin < 4,
  NODELAYHIST=0;
end

if nargin < 5,
  SQRTTIME=1;
end

if nargin < 6,
  PLOTPOL2=0;
else
  PLOTPOL2=1;
end

if PLOTPOL2,
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
linecols = {'r', 'g'};
FONTSIZE=9;
LINEWIDTH=2;

titles = {'pre-mRNA (input)', 'mRNA'};

delayI = 5;
settings = gpnddisimExtractParamTransformSettings(m);

t_len = m.t{1}(end) - m.t{1}(1);

switch SQRTTIME,
  case 2,
    t_pred = linspace(m.t{1}(1), m.t{1}(end), 101)';
    t_tick = m.t{1} - min(m.t{1});
    t_plot = t_pred - min(t_pred);
    t_tickplot = t_tick;
  case 1,
    t_pred = (((0:100)/100*sqrt(t_len)).^2 + 300)';
    t_tick = m.t{1} - min(m.t{1});
    t_plot = sqrt(t_pred - min(t_pred));
    t_tickplot = sqrt(t_tick);
  case 0,
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
    if PLOTPOL2,
      subplot(PLOTROWS, PLOTCOLS, k+1);
    else
      subplot(PLOTROWS, PLOTCOLS, k);
    end
  else
    if PLOTPOL2,
      subplot(PLOTROWS, PLOTCOLS, 3*k+1:3*k+2);
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
      maxbound = ceil(1.2 * max(m.y(J)));
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
    plot(t_tickplot, m.y(J), 'bo')
   case 2,
    errorbar(t_tickplot, m.y(J), 2*sqrt(m.kern.comp{2}.comp{2}.fixedvariance), 'bo')
  end
  hold off
  axis([0, max(t_plot), minbound, maxbound]);
  set(gca, 'XTick', t_tickplot);
  set(gca, 'XTickLabel', t_tick)
  set(gca, 'FontSize', FONTSIZE);
  ylabel(titles{k});
  if k==2,
    xlabel('t (min)');
  end
end

if PLOTPOL2,
  if NODELAYHIST,
    subplot(PLOTROWS, PLOTCOLS, 1);
  else
    subplot(PLOTROWS, PLOTCOLS, 1:2);
  end
  plot(t_tickplot, pol2data, 'bo-')
  axis([0, max(t_plot), 0, 1.1*(max(pol2data)+1e-3)])
  set(gca, 'XTick', t_tickplot);
  set(gca, 'XTickLabel', t_tick)
  set(gca, 'FontSize', FONTSIZE);
  ylabel('Pol II');
end

if ~NODELAYHIST,
  subplot(PLOTROWS, PLOTCOLS, PLOTCOLS*(1:PLOTROWS));
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
    subplot(PLOTROWS, PLOTCOLS, 1);
    title(gene_name);
    subplot(PLOTROWS, PLOTCOLS, PLOTROWS);
  else
    subplot(PLOTROWS, PLOTCOLS, 1:2);
    title(gene_name);
    subplot(PLOTROWS, PLOTCOLS, PLOTCOLS*(1:PLOTROWS));
  end
end
