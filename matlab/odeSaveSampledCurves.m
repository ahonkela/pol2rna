% analyse_ode_full

ID = '2015-05-21_spl1';

load(sprintf('ode_mcmc_%s_summary.mat', ID), 'res');

if ~exist('myI', 'var'),
  error('myI undefined');
end

pol2rnaseqToolboxes;

myI = intersect(myI, 1:size(res, 1));

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

curves = cell(1, size(res, 1));
resvariances = zeros(1, size(res, 1));
datavariances = zeros(1, size(res, 1));

fname = sprintf('ode_mcmc_%s_curves_%d.mat', ID, min(myI));

for k=myI,
  fprintf('Doing gene %d/%d...\n', find(k==myI), length(myI));
  geneI = find(strcmp(res{k,1}.gene_name, d.gene_name));
  mysamples = zeros(400, 6);
  for l=1:4,
    mysamples((1:100)+(l-1)*100, :) = res{k,l}.samples(101:end, :);
  end
  yout = zeros(length(t_plot), size(mysamples, 1));
  mydataVals1 = d.dataVals1(:, geneI) - min(d.dataVals1(:, geneI));
  for n=1:size(mysamples, 1),
    params = num2cell(odeTransformParams(mysamples(n,:)));
    yout(:, n) = odeSimulate(mydataVals1, t_tick, t_pred, params, 0);
  end
  mydataVals2 = d.dataVals2(:, geneI);
  rnascale = max(mydataVals2) / 10;
  mydataVals2 = mydataVals2 / rnascale;

  resvariances(k) = mean(mean((yout(t_tickplot, :) - repmat(mydataVals2, 1, 400)).^2));
  datavariances(k) = var(mydataVals2);
  curves{k} = prctile(yout, [2.5, 50, 97.5], 2);
  safeSave(fname, curves, resvariances, datavariances);
end
