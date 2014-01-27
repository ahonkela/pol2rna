%mypath = '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-08';
%mypath = '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-14';
mypaths = {'~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-08',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-08',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-14',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-14',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-14',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-14',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-17',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-17',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-22',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-22',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-22',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-22',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-23',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-23',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-23',
	  '~/Dropbox/projects/pol2rnaseq/hmc_synthetic_results_2013-05-23'};
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
  fnames{k,l,1} = sprintf('Synthetic_%d_samples_synth_2013-05-07_init%%d.mat', i);
  fnames{k,l,2} = sprintf('Synthetic_%d_samples_synth_nodelay_2013-05-07_init%%d.mat', i);
  fnames{k,l,3} = sprintf('Synthetic_%02d_samples_synth_2_2013-05-14_init%%d.mat', i);
  fnames{k,l,4} = sprintf('Synthetic_%02d_samples_synth_2_nodelay_2013-05-14_init%%d.mat', i);
  fnames{k,l,5} = sprintf('Synthetic_%02d_samples_synth_3_2013-05-14_init%%d.mat', i);
  fnames{k,l,6} = sprintf('Synthetic_%02d_samples_synth_3_nodelay_2013-05-14_init%%d.mat', i);
  fnames{k,l,7} = sprintf('Synthetic_%02d_samples_synth_1_2013-05-17_init%%d.mat', i);
  fnames{k,l,8} = sprintf('Synthetic_%02d_samples_synth_1_nodelay_2013-05-17_init%%d.mat', i);
  fnames{k,l,9} = sprintf('Synthetic_%02d_samples_synth_1_2013-05-22_init%%d.mat', i);
  fnames{k,l,10} = sprintf('Synthetic_%02d_samples_synth_1_nodelay_2013-05-22_init%%d.mat', i);
  fnames{k,l,11} = sprintf('Synthetic_%02d_samples_synth_pol2fix_1_2013-05-22_init%%d.mat', i);
  fnames{k,l,12} = sprintf('Synthetic_%02d_samples_synth_pol2fix_1_nodelay_2013-05-22_init%%d.mat', i);
  fnames{k,l,13} = sprintf('Synthetic_%02d_samples_synth_1_2013-05-23_init%%d.mat', i);
  fnames{k,l,14} = sprintf('Synthetic_%02d_samples_synth_1_nodelay_2013-05-23_init%%d.mat', i);
  fnames{k,l,15} = sprintf('Synthetic_%02d_samples_synth_pol2fix_1_2013-05-23_init%%d.mat', i);
  fnames{k,l,16} = sprintf('Synthetic_%02d_samples_synth_pol2fix_1_nodelay_2013-05-23_init%%d.mat', i);
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
delaymeans = cellfun(@(x) x(5), parammeans(:, :, [1 3 5 7 9 11 13 15]));
settings = gpnddisimExtractParamTransformSettings(r.m);

decays = sigmoidabTransform(decaymeans, 'atox', settings{3});
delays = sigmoidabTransform(delaymeans, 'atox', settings{5});
