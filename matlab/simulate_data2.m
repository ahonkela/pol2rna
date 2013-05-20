% Generate synthetic data from the model

r = load('~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-08/Synthetic_1_samples_synth_2013-05-07_init1.mat');
load('simulated_data.mat')

randn('state', 42);
t_gen{4} = t_gen{1};
N_GEN = 24;
par_gen = cell(N_GEN, 1);
newrna = zeros(N_GEN, 1, length(t));
newpol2 = zeros(N_GEN, 1, length(t));
for k=1:N_GEN,
  par_gen{k} = 2*randn(1, 10);
  m2 = modelExpandParam(r.m, p);
  obs = mvnrnd([ones(10, 1); zeros(10, 1)], m2.K / (0.1*max(max(m2.K))));
  obs(1:10) = obs(1:10) / sqrt(mean(obs(1:10).^2));
  obs(11:20) = obs(11:20) / sqrt(mean(obs(11:20).^2));
  newrna(k,1,:) = obs(1:10);
  newpol2(k,1,:) = obs(11:20);
end

pol2data{4} = newpol2;
rnadata{4} = newrna;

save simulated_data.mat pol2data rnadata Dvals Deltavals t_gen par_gen
