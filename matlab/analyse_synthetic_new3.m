mypaths = {'~/mlprojects/pol2rnaseq/matlab/hmc_results/joint/', ...
           '~/mlprojects/pol2rnaseq/matlab/hmc_results/joint_nodelay/',...
           '~/mlprojects/pol2rnaseq/matlab/hmc_results/joint/', ...
           '~/mlprojects/pol2rnaseq/matlab/hmc_results/joint_nodelay/',...
          };
sdata = load('simulated_data_2013-08-09.mat');

NINITS = 4;
variables = {[1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:9], [1:4, 6, 8:9]};

samples = cell(6, 4, length(mypaths));
fnames = cell(6, 4, length(mypaths));
epsilons = zeros(6, 4, length(mypaths));

datasize = size(sdata.rnadata{1});
for i=1:prod(datasize(1:2)),
  [k, l] = ind2sub(datasize(1:2), i);
  fnames{k,l,1} = sprintf('Synthetic_%02d_samples_synth_new_2013-08-30a_data1_unif0.mat', i);
  fnames{k,l,2} = sprintf('Synthetic_%02d_samples_synth_new_nodelay_2013-08-30a_data1_unif0.mat', i);
  fnames{k,l,3} = sprintf('Synthetic_%02d_samples_synth_new_2013-08-30a_data2_unif0.mat', i);
  fnames{k,l,4} = sprintf('Synthetic_%02d_samples_synth_new_nodelay_2013-08-30a_data2_unif0.mat', i);
end

%d = dir([mypath '/*.mat']);
%filenames = {};
%[filenames{1:length(d),1}] = deal(d.name);

for k=1:length(sdata.Dvals),
  for l=1:length(sdata.Deltavals),
    for n=1:size(samples, 3),
      samples{k,l,n} = zeros(500, length(variables{n}), NINITS);
      r = load([mypaths{n}, '/', fnames{k,l,n}]);
      epsilons(k, l, n) = r.options.epsilon;
      if r.samples_done > 0,
        for m=1:length(r.HMCsamples),
          samples{k,l,n}(:, :, m) = r.HMCsamples{m}(501:end,variables{n});
        end
      end
    end
  end
end

N = size(samples{1,1,1}, 1);
M = size(samples{1,1,1}, 3);

means = cellfun(@(x) squeeze(mean(x, 1)), samples, 'UniformOutput', false);
stds = cellfun(@(x) squeeze(std(x, 0, 1)), samples, 'UniformOutput', false);
W = cellfun(@(x) mean(x.^2, 2), stds, 'UniformOutput', false);
B = cellfun(@(x) N*var(x, 0, 2), means, 'UniformOutput', false);
varHatPlus = cellfun(@(W, B) (N-1)/N * W + 1/N * B, W, B, 'UniformOutput', false);
Rhat = cellfun(@(a, b) sqrt(a ./ b), varHatPlus, W, 'UniformOutput', false);

difficult = cellfun(@(x) any(x > 1.2), Rhat);

parammeans = cellfun(@(x) squeeze(mean(x, 2)), means, 'UniformOutput', false);
decaymeans = cellfun(@(x) x(3), parammeans);
sensmeans = cellfun(@(x) x(4), parammeans);
delaymeans = cellfun(@(x) x(5), parammeans(:, :, 1:2:end));
settings = gpnddisimExtractParamTransformSettings(r.m);

decays = sigmoidabTransform(decaymeans, 'atox', settings{3});
delays = sigmoidabTransform(delaymeans, 'atox', settings{5});
