seeds = 11:14;
id = '2015-05-15';
MYP = [2.5 25 50 75 97.5];
NPARAMS = 6;

d = dir('ode_mcmc_results/*.mat');
genes0 = cell(size(d));
for k=1:length(d),
  t = strsplit(d(k).name, '_');
  genes0{k} = t{1};
end
genes = unique(genes0);

res = cell(length(genes), length(seeds));
means = zeros(length(genes), length(seeds), NPARAMS);
medians = zeros(length(genes), length(seeds), NPARAMS);
stds = zeros(length(genes), length(seeds), NPARAMS);
genemedians = zeros(length(genes), NPARAMS);
genestds = zeros(length(genes), NPARAMS);
geneprcts = zeros(length(genes), NPARAMS, length(MYP));
genells = zeros(length(genes), 1);
for k=1:length(genes),
  if mod(k, 10) == 0
    fprintf('Doing gene %d/%d\n', k, length(genes));
  end
  mysamples = zeros(length(seeds)*100, NPARAMS);
  mylls = zeros(length(seeds)*100, 1);
  for l=1:length(seeds),
    fname = sprintf('ode_mcmc_results/%s_samples_%s_seed%d.mat', genes{k}, id, seeds(l));
    res{k,l} = load(fname);
    means(k, l, :) = mean(res{k,l}.samples(101:end, :));
    medians(k, l, :) = median(res{k,l}.samples(101:end, :));
    stds(k, l, :) = std(res{k,l}.samples(101:end, :));
    mysamples((1:100)+(l-1)*100, :) = res{k,l}.samples(101:end, :);
    mylls((1:100)+(l-1)*100) = res{k,l}.ll(101:end);
  end
  genemedians(k, :) = median(mysamples);
  genestds(k, :) = std(mysamples);
  geneprcts(k, :, :) = prctile(mysamples, MYP)';
  genells(k) = mean(mylls);
end

N = size(res{k,l}.samples, 1) - 100;
W = squeeze(mean(stds.^2, 2));
B = squeeze(var(means, [], 2));
varHatPlus = (N-1)/N * W + 1/N * B;
Rhat = sqrt(varHatPlus./W);

d = dir(sprintf('ode_mcmc_summaries/ode_mcmc_%s_curves*.mat', id));
fnames = cellfun(@(x) ['ode_mcmc_summaries/' x], {d.name}, 'UniformOutput', 0);
r = merge_files(fnames);

truemedians = odeTransformParams(genemedians);
trueprcts = zeros(size(geneprcts));
for l=1:length(MYP),
  trueprcts(:, :, l) = odeTransformParams(geneprcts(:, :, l));
end
fp = fopen(sprintf('results/ctd_delays_%s.txt', id), 'w');
fprintf(fp, 'gene\tctd_2.5%%\tctd_25%%\tctd_50%%\tctd_75%%\tctd_97.5%%\tavell\tctd_iqr\tctd_resvarfrac\n');
for k=1:length(genes),
  fprintf(fp, '%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', genes{k}, squeeze(trueprcts(k, 4, :)), genells(k), abs(diff(trueprcts(k, 4, [4,2]))), r.resvariances(k)/r.datavariances(k));
end
fclose(fp);



% genes = { 'ENSG00000163634',
% 'ENSG00000115641',
% 'ENSG00000197451',
% 'ENSG00000185090',
% 'ENSG00000126756',
% 'ENSG00000171067',
% 'ENSG00000166166',
% 'ENSG00000160058',
% 'ENSG00000179913',
% 'ENSG00000197375',
% 'ENSG00000214530',
% 'ENSG00000197170',
% 'ENSG00000107815',
% 'ENSG00000115875',
% 'ENSG00000132661',
% 'ENSG00000110092',
% 'ENSG00000187145',
% 'ENSG00000196139',
% 'ENSG00000173890',
% 'ENSG00000183578',
% 'ENSG00000064195',
% 'ENSG00000181610',
% 'ENSG00000170442',
% 'ENSG00000169255',
% 'ENSG00000126524',
% 'ENSG00000258289',
% 'ENSG00000167173',
% 'ENSG00000131143',
% 'ENSG00000092010',
% 'ENSG00000137038',
% 'ENSG00000120533',
% 'ENSG00000181026',
% 'ENSG00000103257',
% 'ENSG00000131051',
% 'ENSG00000062725',
% 'ENSG00000167767',
% 'ENSG00000072210',
% 'ENSG00000128567',
% 'ENSG00000168140',
% 'ENSG00000198910',
% 'ENSG00000099992',
% 'ENSG00000056736',
% 'ENSG00000125691',
% 'ENSG00000142864',
% 'ENSG00000164125',
% 'ENSG00000048162',
% 'ENSG00000023445',
% 'ENSG00000123066',
% 'ENSG00000173227',
% 'ENSG00000179388',
% 'ENSG00000136603',
% 'ENSG00000164684',
% 'ENSG00000143126',
% 'ENSG00000173821',
% 'ENSG00000164938',
% 'ENSG00000109452',
% 'ENSG00000182022',
% 'ENSG00000165272',
% 'ENSG00000121671',
% 'ENSG00000154767',
% 'ENSG00000103061',
% 'ENSG00000197622',
% 'ENSG00000112701',
% 'ENSG00000120688',
% 'ENSG00000166949',
% 'ENSG00000119285',
% 'ENSG00000254635',
% 'ENSG00000143553',
% 'ENSG00000092969',
% 'ENSG00000127084',
% 'ENSG00000006652',
% 'ENSG00000099622',
% 'ENSG00000144228',
% 'ENSG00000156127',
% 'ENSG00000135052',
% 'ENSG00000248099',
% 'ENSG00000073282',
% 'ENSG00000166483',
% 'ENSG00000071242',
% 'ENSG00000113971',
% 'ENSG00000119318',
% 'ENSG00000096968',
% 'ENSG00000240498',
% 'ENSG00000106003',
% 'ENSG00000162949',
% 'ENSG00000124831',
% 'ENSG00000177119',
% 'ENSG00000168209',
% 'ENSG00000258890',
% 'ENSG00000140548',
% 'ENSG00000163993',
% 'ENSG00000138119',
% 'ENSG00000173848',
% 'ENSG00000135245',
% 'ENSG00000103035',
% 'ENSG00000109321',
% 'ENSG00000184012',
% 'ENSG00000157483',
% 'ENSG00000166888',
% 'ENSG00000143379' };
