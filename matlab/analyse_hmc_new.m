resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/joint/';
id = '2013-08-30';
%id = '2013-11-05';
%resultdir = '/share/synergy/analyses/hmc_results/joint/';

d = dir([resultdir '*_samples_' id '_unif0.mat']);

filenames = {};
[filenames{1:length(d),1}] = deal(d.name);
%filenames = filenames(cellfun('length', filenames) == 44);

find_good_pol2;
[~, A, B] = intersect(pol2ensg, cellfun(@(x) x(1:15), filenames, 'UniformOutput', false));
filenames = filenames(B);

%N_GOOD = 5;

Isampl = 501:1000;
Isampl_thin = 501:10:1000;

MYP = [2.5 25 50 75 97.5];

numgenes = length(filenames);

%parI = [1:4, 6, 8:10];
parI = [1:6, 8:10];
mygene = '';
genes = cell(length(filenames), 1);
means = zeros(length(filenames), length(parI));
stds = zeros(length(filenames), length(parI));
medians = zeros(length(filenames), length(parI));
thetameans = zeros(length(filenames), 1);
prcts = zeros(length(filenames), length(parI), length(MYP));
curI = 1:length(Isampl_thin);
Iincr = length(curI);
h = cell(length(filenames), 1);
for k=1:length(filenames),
  fprintf('%d/%d\n', k, length(filenames));
  r = load([resultdir, filenames{k}]);
  genes{k} = r.gene_name;
  if r.finished,
    hk = cat(3, r.HMCsamples{:});
    hk = hk(Isampl, parI, :);
    hk = permute(hk, [1 3 2]);
    sz = size(hk);
    hk = reshape(hk, [sz(1)*sz(2), sz(3)]);
    prcts(k, :, :) = prctile(hk, MYP)';
    means(k, :) = mean(hk);
    medians(k, :) = median(hk);
    stds(k, :) = std(hk);
  end
end

% p30 = zeros(length(filenames), 1);
% bound30 = sigmoidabTransform(30, 'xtoa', [0, 299]);

% for k=1:length(h),
%   p30(k) = mean(h{k}(:, 5) > bound30);
% end

% baddata = zeros(numgenes, 1);
% Rhat = zeros(numgenes, length(parI));
% thetatruemeans = zeros(numgenes, 1);
% for k=1:numgenes,
%   I = (k-1)*N_GOOD + (1:N_GOOD);
%   S = sum(bsxfun(@minus, means(I,:), median(means(I,:))) .^ 2, 2);
%   J = I(S < 100 & sum(stds(I, :), 2) > 0.1);
%   if length(J) < N_GOOD,
%     baddata(k) = 1;
%     disp(J)
%   end
%   N = length(Isampl);
%   M = length(I);
%   W = mean(stds(I,:).^2);
%   B = N*var(means(I,:));
%   varHatPlus = (N-1)/N * W + 1/N * B;
%   Rhat(k,:) = sqrt(varHatPlus./W);
%   thetatruemeans(k) = mean(thetameans(I));
% end

% K = max(Rhat, [], 2) < 1.2;
% goodgenes = genes(find(K)*5);
% fp = fopen('finished_genes_latest.txt', 'w');
% fprintf(fp, '%s\n', goodgenes{:});
% fclose(fp);

% Kbad = ~(max(Rhat, [], 2) < 1.2);
% badgenes = genes(find(Kbad)*5);
% fp = fopen('unfinished_genes_latest.txt', 'w');
% fprintf(fp, '%s\n', badgenes{:});
% fclose(fp);

% save analysis_dump

% for k=find(max(Rhat, [], 2) > 1.2)',
%   I = (k-1)*5 + (1:5);
%   fprintf('k=%d\n', k);
%   J = Rhat(k,:) > 1.2;
%   disp(Rhat(k,:));
%   disp(means(I,:));
%   disp(stds(I,:));
%   % for k=1:length(I),
%   %   subplot(5, 1, k);
%   %   plot(h{I(k)}(:,J));
%   % end
%   % pause
% end

save(['results/result_summary_' id], 'means', 'stds', 'prcts', 'genes');

pr = load(['results/early_profiles_' id '.mat']);
diff1 = max(pr.mu(:, 1:20), [], 2) - min(pr.mu(:, 1:20), [], 2);
diff2 = max(pr.mu(:, 21:30), [], 2) - min(pr.mu(:, 21:30), [], 2);
qrtl = prctile(diff1-diff2, [25, 75]);
bound = qrtl(2) + 1.5 * (qrtl(2) - qrtl(1));
early_dev = (diff1 - diff2 - qrtl(2)) / (qrtl(2) - qrtl(1));

diff1 = (max(pr.mu(:, 1:20), [], 2) - ...
         min(pr.mu(:, 1:20), [], 2)) ./ max(pr.mu, [], 2);
diff2 = (max(pr.mu(:, 21:30), [], 2) - ...
         min(pr.mu(:, 21:30), [], 2)) ./ max(pr.mu, [], 2);
begdev10b = diff1 - diff2;

[~, A, B] = intersect(genes, pr.genes);
assert(all(strcmp(genes(A), pr.genes(B))));
early_devs = zeros(length(genes), 1);
early_devs(A) = early_dev(B);
begdevs = zeros(length(genes), 1);
begdevs(A) = begdev10b(B);

delaypcts = sigmoidabTransform(squeeze(prcts(:,5,:)), 'atox', [0, 299]);
fp = fopen(['results/hmc_results_to_browser_' id '.txt'], 'w');
fprintf(fp, 'gene 2.5%% 25%% 50%% 75%% 97.5%% begdev\n');
for k=1:length(genes),
  if ~all(delaypcts(k,:) == 149.5),
    fprintf(fp, '%s %f %f %f %f %f %f\n', genes{k}, delaypcts(k,:), begdevs(k));
  end
end
fclose(fp);
