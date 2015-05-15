files = dir('ode_mcmc_results/*.mat');
fnames = cell(size(files));
for k=1:length(fnames),
  fnames{k} = ['ode_mcmc_results/' files(k).name];
end

I = intersect(myI, 1:length(fnames));

odeAddLikelihoods(fnames(I));
