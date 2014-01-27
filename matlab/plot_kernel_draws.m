% Generate synthetic data from the model

FONTSIZE=9;

r = load('~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-08/Synthetic_5_samples_synth_nodelay_2013-05-07_init1.mat');
load('simulated_data.mat')

t = (1:2:200)' + 300;

randn('state', 44);
t_gen{4} = t_gen{1};
N_GEN = 3;
par_gen = cell(N_GEN, 1);
newrna = zeros(N_GEN, length(t));
newpol2 = zeros(N_GEN, length(t));
m0 = gpnddisimConvertTransformSettings(r.m);
k = 1;
%par_gen{k} = 2*randn(1, 10);
%par_gen{k}(5) = -100;
%m2 = modelExpandParam(m0, par_gen{k});
  %obs = mvnrnd([ones(10, 1); zeros(10, 1)], m2.K / (0.1*max(max(m2.K))));

decays = [-4, -2, 0];
filenames = {'kern_samples_1.eps', 'kern_samples_2.eps', 'kern_samples_3.eps'};

for decI=1:length(decays),
  par = gpnddisimExtractParam(m0);
  % inverse width
  par(1) = 1;
  % variance
  par(2) = -6;
  % decay
  par(3) = decays(decI);
  % sensitivity
  par(4) = 3;
  %par(4) = -5  + 2*decays(decI);
  m0 = gpnddisimExpandParam(m0, par);
  kern = m0.kern.comp{1};
  K = kernCompute(kern, {t, t});
  Kstore{decI} = K;
  
  randn('state', 43);
  
  for k=1:N_GEN,
    means = zeros(size(K(:, 1)));
    %means(1:length(t)) = randn(1);
    %means(length(t)+1:end) = 
    obs = mvnrnd(means, K);
    newpol2(k,:) = obs(1:length(t));
    newrna(k,:) = obs(end-length(t)+1:end);
  end
  
  figure(decI);
  subplot(2, 1, 1)
  plot(squeeze(newpol2)')
  set(gca, 'XTick', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'FontSize', FONTSIZE);
  ylabel('Pol II density')
  subplot(2, 1, 2)
  plot(squeeze(newrna)')
  set(gca, 'XTick', []);
  set(gca, 'YTickLabel', []);
  set(gca, 'FontSize', FONTSIZE);
  ylabel('mRNA abundance')
  set(gcf, 'PaperPosition', [0 0 4 3]);
  %print('-depsc2', filenames{decI});
end
