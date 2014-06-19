SAMPLEDIR='joint/';
FILESTEM='_hmc_pol2';
FORMAT='eps';
SQRTTIME=0;
FILESPEC='_samples_2013-08-30_unif0.mat';
PLOTPREMRNA=0;
PREMRNAMODEL=0;
SHORTHIST=1;

genes = {'ENSG00000183484',
         'ENSG00000198807',
         'ENSG00000133816',
         'ENSG00000173227',
         'ENSG00000019549',
         'ENSG00000064195',
         'ENSG00000170027',
         'ENSG00000139116'};

for k=1:length(genes),
  fprintf('Pass 1 %d/%d\n', k, length(genes));
  plot_sampled_predictions_file(genes{k}, SAMPLEDIR, FILESTEM, FORMAT, ...
                                SQRTTIME, FILESPEC, PLOTPREMRNA, ...
                                PREMRNAMODEL, SHORTHIST);
end

SAMPLEDIR='joint/';
FILESTEM='_hmc_premrna';
FORMAT='eps';
SQRTTIME=0;
FILESPEC='_samples_2013-11-05_unif0.mat';
PLOTPREMRNA=0;
PREMRNAMODEL=1;
SHORTHIST=1;

for k=1:length(genes),
  fprintf('Pass 2 %d/%d\n', k, length(genes));
  plot_sampled_predictions_file(genes{k}, SAMPLEDIR, FILESTEM, FORMAT, ...
                                SQRTTIME, FILESPEC, PLOTPREMRNA, ...
                                PREMRNAMODEL, SHORTHIST);
end
