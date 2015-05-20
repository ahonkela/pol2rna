% analyse_ode_full

FONTSIZE=6;
LINEWIDTH=1.5;

d = load('ode_data_2015-05-07.mat');

t_pred = [linspace(0, 157.5, 64), exp(linspace(log(160), log(1280), 61))]';
t_tick = d.timevector - min(d.timevector);
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

curves = cell(1, size(res, 2));
resvariances = zeros(1, size(res, 2));
datavariances = zeros(1, size(res, 2));

for k=1:4,
  geneI = find(strcmp(res{k,1}.gene_name, d.gene_name));
  mysamples = zeros(400, 6);
  for l=1:4,
    mysamples((1:100)+(l-1)*100, :) = res{k,l}.samples(101:end, :);
  end
  yout = zeros(length(t_plot), size(mysamples, 1));
  for n=1:size(mysamples, 1),
    params = num2cell(odeTransformParams(mysamples(n,:)));
    yout(:, n) = odeSimulate(d.dataVals1(:,geneI), ...
                             t_tick, t_pred, ...
                             params, 0);
  end
  mydataVals2 = d.dataVals2(:, geneI);
  rnascale = max(mydataVals2) / 10;
  mydataVals2 = mydataVals2 / rnascale;

  resvariances(k) = mean(mean((yout(t_tickplot, :) - repmat(mydataVals2, 1, 400)).^2));
  datavariances(k) = var(mydataVals2);
  curves{k} = prctile(yout, [2.5, 50, 97.5], 2);
end





for k=1:4,
  geneI = find(strcmp(res{k,1}.gene_name, d.gene_name));
  mydataVals2 = d.dataVals2(:, geneI);
  myrnaVars = d.rnaVars(:, geneI);
  rnascale = max(mydataVals2) / 10;
  mydataVals2 = mydataVals2 / rnascale;
  myrnaVars = myrnaVars / (rnascale.^2);

  minbound = 0;
  maxbound = ceil(1.1 * max(mydataVals2 + 2*sqrt(myrnaVars)));

  mycurve = curves{k};
  mycurve(mycurve(:, 1) < minbound, 1) = minbound;
  mycurve(mycurve(:, 3) > maxbound, 3) = maxbound;
  mycol = [0.9 1.0 0.9];

  figure(k);
  h=fill([t_plot; t_plot(end:-1:1)], [mycurve(:,1); flipud(mycurve(:,3))],mycol);
  set(h,'EdgeColor',mycol);
  hold on
  plot(t_plot, mycurve(:,2), 'g', 'LineWidth', LINEWIDTH);
  errorbar(t_tickplot, mydataVals2, ...
           2*sqrt(myrnaVars), 'bx')
  hold off
  axis([0, max(t_plot), minbound, maxbound]);
  set(gca, 'XTick', t_tickplot);
  set(gca, 'XTickLabel', t_ticklabels)
  set(gca, 'FontSize', FONTSIZE);
  v = axis;
  text(v(1)-0.22*(v(2)-v(1)), mean(v(3:4)), 'mRNA', 'HorizontalAlignment','center', 'Rotation', 90, 'FontSize', FONTSIZE)
  v = axis;
  text(mean(v(1:2)), v(3)-0.26*(v(4)-v(3)), 't (min)', 'HorizontalAlignment','center', 'FontSize', FONTSIZE);
end
