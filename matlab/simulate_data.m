datapath = '~/projects/pol2rnaseq/data/';

TIPARP='ENSG00000163659';
load([datapath, 'pol2_summaryseries_2013_01_02.mat']);
load([datapath, 'bininfo_dec2012_corrected.mat'], 'bininfo');

TIPIND = find(bininfo(:, 5) == ensg2int(TIPARP));
t = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280];
t_test = 0:2.5:1280;

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

Dvals = [0.32, 0.16, 0.08, 0.04, 0.02, 0.01];
funcs = zeros(length(Dvals), length(t_test));

for k=1:length(Dvals),
  S = 0.01;
  B = 0.001;
  D = Dvals(k);
  Delta = 0;
  m0 = 0.001;
  t0 = 0;

  m = m0 * exp(D*(t0-t_test)) + B/D * (1-exp(D*(t0-t_test))) + ...
      S*exp(-D*t_test) .* cumtrapz(exp(D*t_test) .* pol2fit');
  funcs(k,:) = m / max(m);
end

plot(sqrt(t_test), funcs);
hold on;
plot(sqrt(t_test), pol2fit / max(pol2fit), 'k--');
plot(sqrt(t), myrna / max(myrna), 'k-.');
hold off
legend('D=0.32', 'D=0.16', 'D=0.08', 'D=0.04', 'D=0.02', 'D=0.01', 'Pol2', 'RNA', 'Location', 'NorthEast')
axis([0 sqrt(1280) 0 1])
set(gca, 'XTick', sqrt(t));
set(gca, 'XTickLabel', t);
%set(gcf, 'PaperPosition', [0 0 14 10]); print -depsc2 tiparp_simulation.eps
