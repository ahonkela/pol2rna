N = 200;
N_all = 30494;
filestem = 'fittingresults_shifted_longerdelay_new_2012-05-07b_%d_200.mat';

results_geneindices = zeros(1, N_all);
results_ensemblids = zeros(1, N_all);
results_loglikelihoods = zeros(N_all, 3);
results_jointmodels = cell(1, N_all);
results_jointtransforminfos = cell(1, N_all);

for k=1:N,
    fprintf('%d/%d\n', k, N);
    res = load(sprintf(filestem, k));
    I = (res.results_geneindices ~= 0);
    results_geneindices(I) = res.results_geneindices(I);
    results_ensemblids(I) = res.results_ensemblids(I);
    results_loglikelihoods(I, :) = res.results_loglikelihoods(I, :);
    results_jointmodels(I) = res.results_jointmodels(I);
    results_jointtransforminfos(I) = res.results_jointtransforminfos(I);
end
