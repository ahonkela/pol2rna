% Generate synthetic data based on TIPARP Pol2 profile
% Antti Honkela, May 2013

datapath = '~/projects/pol2rnaseq/data/';

TIPARP='ENSG00000163659';
load([datapath, 'pol2_summaryseries_2013_01_02.mat']);
load([datapath, 'bininfo_dec2012_corrected.mat'], 'bininfo');

TIPIND = find(bininfo(:, 5) == ensg2int(TIPARP));
t = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280];
t_test = 0:2.5:1280;
t_gen = {[0, 5, 10, 20, 40, 80, 160, 320, 640, 1280], ...
        [0 5 10:10:240 280:40:640], ...
        [0, 10, 20, 40, 60, 80, 120, 160, 240, 320]};

if ~exist('r', 'var'),
  r = load([datapath, 'info_gene_mean_var.mat']);
end

myrna = r.mu(find(strcmp(r.geneID, TIPARP)), :);

mypol = pol2_summaryseries(TIPIND, :);

% Fit GP to pol2
if ~exist('pol2fit', 'var'),
  optionsGP = gpOptions('ftc');
  
  % The kernel used is a sum of an rbf, a bias and a white term. Bias and
  % white help numerical stability during optimisation (this will also
  % result in bigger error-bars in the end). Also try 'lin' or 'matern32'
  % or even 'rbfperiodic' for the first component of the kernel.
  optionsGP.kern = {'mlp','white'};
  %optionsGP.scale2var1 = true;     % Scale outputs to variance 1.
  modelGP = gpCreate(1, 1, sqrt(t'), mypol', optionsGP);
  modelGP = gpOptimise(modelGP, 1, 1000);
  % The model is now optimised. The following function just does standard
  % GP prediction based on the optimised model and the test inputs.
  [pol2fit, yVar] = gpPosteriorMeanVar(modelGP, sqrt(t_test'));
  pol2fit = pol2fit - min(pol2fit);
end
% end GP fit

Dvals = log(2) ./ [2 4 8 16 32 64]; %[0.32, 0.16, 0.08, 0.04, 0.02, 0.01];
Deltavals = [0 4 8 12];
funcs = zeros(length(Dvals), length(Deltavals), length(t_test));

for k=1:length(Dvals),
  for l=1:length(Deltavals),
    D = Dvals(k);
    S = 0.03;
    B = 0.005;
    Delta = Deltavals(l);
    m0 = 0.008 / D;
    t0 = 0;

    pol2shift = [pol2fit(1)*ones(Delta, 1); pol2fit(1:end-Delta)];
  
    m = m0 * exp(D*(t0-t_test)) + B/D * (1-exp(D*(t0-t_test))) + ...
	S*exp(-D*t_test) .* cumtrapz(exp(D*t_test) .* pol2shift');
    funcs(k,l,:) = 2*m / sqrt(mean(m.^2));
  end
end

t_indices = zeros(size(t));
for k=1:length(t),
  t_indices(k) = find(t_test == t(k));
end

plot(sqrt(t_test(t_indices)), squeeze(funcs(:, 1, t_indices)));
hold on;
plot(sqrt(t_test), pol2fit / max(pol2fit), 'k--');
plot(sqrt(t), myrna / max(myrna), 'k-.');
hold off
legend('t1/2=2', 't1/2=4', 't1/2=8', 't1/2=16', 't1/2=32', 't1/2=64', 'Pol2', 'RNA', 'Location', 'NorthEast')
axis([0 sqrt(1280) 0 5])
set(gca, 'XTick', sqrt(t));
set(gca, 'XTickLabel', t);
%set(gcf, 'PaperPosition', [0 0 14 10]); print -depsc2 tiparp_simulation.eps

rnadata = cell(size(t_gen));
pol2data = cell(size(t_gen));
for k=1:length(t_gen),
  t_indices = zeros(size(t_gen{k}));
  for l=1:length(t_indices),
    t_indices(l) = find(t_test == t_gen{k}(l));
  end

  pol2data{k} = pol2fit(t_indices);
  rnadata{k} = funcs(:, :, t_indices);
  randn('state', 42);
  pol2noise = 0.1 * randn(size(pol2data{k}));
  rnanoise = 0.1 * randn(1, 1, length(t_indices));
  pol2data{k} = pol2data{k} + pol2noise;
  rnadata{k} = rnadata{k} + repmat(rnanoise, [size(rnadata{k}, 1), size(rnadata{k}, 2), 1]);
end

save simulated_data_2013-08-09.mat pol2data rnadata Dvals Deltavals t_gen
