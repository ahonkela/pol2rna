[a, names] = gpnddisimExtractParamTransformSettings(m);
[foo, names] = gpnddisimExtractParam(m);
b = [a(1:6), a(8:10)];
names = [names(1:6), names(8:10)];
for k=1:length(names),
  I = strfind(names{k}, ',');
  if ~isempty(I),
    names{k} = names{k}(1:I-1);
  end
end
for k=1:6,
  names{k} = names{k}(9:end);
end

if 0,
for i = 1:length(genes),
  fprintf('%d/%d\n', i, length(genes));
  clf;
  plot_hists(h{i}, b, names, 20)
  subplot(3, 3, 8);
  v = axis;
  text(mean(v(1:2)), v(3) - 0.25*v(4), genes{i}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold')
  print('-dpng', '-r100', sprintf('plots/%s.png', genes{i}));
end
end

numbers = zeros(size(genes));
indices = zeros(size(genes));
for k=1:length(genes),
  numbers(k) = str2double(genes{k}(5:end));
  indices(k) = find(bininfo(:, 5) == numbers(k));
end

if 0,
for k = 1:length(genes),
  fprintf('%d/%d\n', k, length(genes));
  clf;
  plot_data(pol_summaryseries(indices(k),:), rna_summaryseries(indices(k),:), genes{k})
  print('-dpng', '-r100', sprintf('data_plots/%s.png', genes{k}));
end
end

means = zeros(size(genes));
medians = zeros(size(genes));
q1 = zeros(size(genes));
q5 = zeros(size(genes));
q10 = zeros(size(genes));
q95 = zeros(size(genes));
for k=1:length(genes),
  fprintf('%d/%d\n', k, length(genes));
  mysamples = sigmoidabTransform(h{k}(:, 5), 'atox', b{5});
  means(k) = mean(mysamples);
  medians(k) = median(mysamples);
  q1(k) = quantile(mysamples, 0.01);
  q5(k) = quantile(mysamples, 0.05);
  q10(k) = quantile(mysamples, 0.10);
  q95(k) = quantile(mysamples, 0.95);
end

fid = fopen('delay_results_2012-01-10.txt', 'w');
fprintf(fid, 'gene mean q01 q05 q10 q50 q95\n');
for k=1:length(genes),
  fprintf(fid, '%s %f %f %f %f %f %f\n', genes{k}, means(k), q1(k), q5(k), q10(k), medians(k), q95(k));
end
fclose(fid);
