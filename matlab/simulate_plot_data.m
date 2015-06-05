% Generate synthetic data based on TIPARP Pol2 profile
% Antti Honkela, May 2013

datapath = '~/projects/pol2rnaseq/data/';

TIPARP='ENSG00000163659';
load([datapath, 'pol2_summaryseries_2013_01_02.mat']);
load([datapath, 'bininfo_dec2012_corrected.mat'], 'bininfo');

TIPIND = find(bininfo(:, 5) == ensg2int(TIPARP));
t = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280];
t_test = 0:2.5:1280;
t_gen = {t_test};

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

Dvals = log(2) ./ [8 32 8]; %[0.32, 0.16, 0.08, 0.04, 0.02, 0.01];
Deltavals = [0 0 8];
funcs = zeros(length(Dvals), length(t_test));

for k=1:length(Dvals),
  D = Dvals(k);
  S = 0.03;
  B = 0.005;
  Delta = Deltavals(k);
  m0 = 0.008 / D;
  t0 = 0;

  params = {B, 1/(2.5^2)*S^2, D, t_test(Delta+1), m0};
  m2 = odeSimulate(pol2fit, t_test, t_test, params);
  %funcs(k,l,:) = 2*m2 / sqrt(mean(m2.^2));
  funcs(k,:) = 4*m2 / max(m2);
end

t_tick = [0, 20, 40, 80, 160, 320];

subplot(2, 1, 1);
plot(t_test, 4*pol2fit / max(pol2fit));
set(gca, 'FontSize', 6);
axis([0 400 0 4])
set(gca, 'YTick', [])
set(gca, 'XTick', t_tick);
set(gca, 'XTickLabel', t_tick);
%xlabel('t (min)')
ylabel('Pol-II activity')
subplot(2, 1, 2);
plot(t_test, funcs);
axis([0 400 0 4])
set(gca, 'FontSize', 6);
set(gca, 'XTick', t_tick);
set(gca, 'XTickLabel', t_tick);
set(gca, 'YTick', []);
xlabel('t (min)')
ylabel('mRNA activity')
set(gcf, 'PaperPositionMode','auto')
set(gcf, 'PaperPosition', [0 0 87/25.4 56/25.4]);
legh= legend('t_{1/2} = 8 min', 't_{1/2} = 32 min', '\Delta = 20 min', 'Location', 'NorthEast');
legh.Position = legh.Position + [0.0001 0.0001 0.0001 0.0001];
print -depsc2 cartoon_simulation.eps
