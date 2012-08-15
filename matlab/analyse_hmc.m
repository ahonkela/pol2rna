resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/joint/';
%resultdir = '/share/synergy/analyses/hmc_results/joint/';

d = dir([resultdir '*.mat']);

filenames = {};
[filenames{1:length(d),1}] = deal(d.name);
% exclude init3, as it seems very unreliable
I = cellfun('isempty', strfind(filenames, 'init3'));
filenames = filenames(I);

ids = zeros(size(filenames));
for k=1:length(filenames),
  ids(k) = str2num(filenames{k}(5:15));
end
goodids = find(accumarray(ids, 1) == 5);
Ifiles = ismember(ids, goodids);

Isampl = 501:1000;
Isampl_thin = 501:10:1000;

MYP = [2.5 25 50 75 9.75];

filenames = filenames(Ifiles);

mygene = '';
genes = cell(length(filenames), 1);
means = zeros(length(filenames), 9);
stds = zeros(length(filenames), 9);
prcts = zeros(length(filenames)/5, 9, length(MYP));
hgene = zeros(length(Isampl)/2, 9);
curI = 1:length(Isampl_thin);
Iincr = length(curI);
%h = cell(length(filenames), 1);
for k=1:length(filenames),
  fprintf('%d/%d\n', k, length(filenames));
  r = load([resultdir, filenames{k}]);
  genes{k} = r.gene_name;
  hk = r.HMCsamples(Isampl, [1:6, 8:10]);
  if ~strcmp(r.gene_name, mygene),
    if k>1,
      prcts(floor(k/5), :, :) = prctile(hgene, MYP)';
    end
    curI = 1:length(Isampl_thin);
    mygene = r.gene_name;
  end
  hgene(curI,:) = r.HMCsamples(Isampl_thin, [1:6, 8:10]);
  curI = curI + Iincr;
  means(k, :) = mean(hk);
  stds(k, :) = std(hk);
end
prcts(end, :, :) = prctile(hgene, MYP)';

% p30 = zeros(length(filenames), 1);
% bound30 = sigmoidabTransform(30, 'xtoa', [0, 299]);

% for k=1:length(h),
%   p30(k) = mean(h{k}(:, 5) > bound30);
% end

baddata = zeros(length(goodids), 1);
Rhat = zeros(length(goodids), 9);
for k=1:length(goodids),
  I = (k-1)*5 + (1:5);
  S = sum(bsxfun(@minus, means(I,:), median(means(I,:))) .^ 2, 2);
  I = I(S < 100 & sum(stds(I, :), 2) > 0.1);
  if length(I) < 5,
    baddata(k) = 1;
    disp(I)
  end
  N = length(Isampl);
  M = length(I);
  W = mean(stds(I,:).^2);
  B = N*var(means(I,:));
  varHatPlus = (N-1)/N * W + 1/N * B;
  Rhat(k,:) = sqrt(varHatPlus./W);
end

K = max(Rhat, [], 2) < 1.2;
goodgenes = genes(find(K)*5);
fp = fopen('finished_genes_latest.txt', 'w');
fprintf(fp, '%s\n', goodgenes{:});
fclose(fp);

save analysis_dump

for k=find(max(Rhat, [], 2) > 1.2)',
  I = (k-1)*5 + (1:5);
  S = sum(bsxfun(@minus, means(I,:), median(means(I,:))) .^ 2, 2);
  I = I(S < 100 & sum(stds(I, :), 2) > 0.1);
  fprintf('k=%d\n', k);
  J = Rhat(k,:) > 1.2;
  disp(Rhat(k,:));
  disp(means(I,:));
  disp(stds(I,:));
  % for k=1:length(I),
  %   subplot(5, 1, k);
  %   plot(h{I(k)}(:,J));
  % end
  pause
end


delaypcts = sigmoidabTransform(squeeze(prcts(K,5,:)), 'atox', [0, 299]);
fp = fopen('delay_prctiles_2012-08-14.txt', 'w');
fprintf(fp, 'gene 2.5%% 25%% 50%% 75%% 97.5%%\n');
for k=1:length(goodgenes),
  fprintf(fp, '%s %f %f %f %f %f\n', goodgenes{k}, delaypcts(k,:));
end
fclose(fp);
