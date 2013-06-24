mypaths = {'~/mlprojects/pol2rnaseq/matlab/hmc_results/joint/',
          '~/mlprojects/pol2rnaseq/matlab/hmc_results/joint_nodelay/'};
sdata = load('simulated_data.mat');

inits = [1,2,4,5, 8];
%inits = [2,4,8];
variables = {[1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10], ...
	     [1:6, 8:10], [1:4, 6, 8:10]};

samples = cell(6, 4, length(mypaths));
fnames = cell(6, 4, length(mypaths));

datasize = size(sdata.rnadata{1});
for i=1:prod(datasize(1:2)),
  [k, l] = ind2sub(datasize(1:2), i);
  fnames{k,l,1} = sprintf('Synthetic_%02d_samples_synth_1_2013-06-21b_init%%d.mat', i);
  fnames{k,l,2} = sprintf('Synthetic_%02d_samples_synth_1_nodelay_2013-06-21b_init%%d.mat', i);
end

%d = dir([mypath '/*.mat']);
%filenames = {};
%[filenames{1:length(d),1}] = deal(d.name);

for k=1:length(sdata.Dvals),
  for l=1:length(sdata.Deltavals),
    for n=1:size(samples, 3),
      samples{k,l,n} = zeros(500, length(variables{n}), length(inits));
      for m=1:length(inits),
        r = load([mypaths{n}, '/', sprintf(fnames{k,l,n}, inits(m))]);
        samples{k,l,n}(:, :, m) = r.HMCsamples(501:end,variables{n});
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
delaymeans = cellfun(@(x) x(5), parammeans(:, :, [1]));
%settings = gpnddisimExtractParamTransformSettings(r.m);

%decays = sigmoidabTransform(decaymeans, 'atox', settings{3});
%delays = sigmoidabTransform(delaymeans, 'atox', settings{5});
