% Plot heatmaps

% -----------------------------------------------------------------
% GP POSTERIORS
% -----------------------------------------------------------------

highlight_genes = {'ENSG00000162949',
		   'ENSG00000019549',   %'ENSG00000168140',
                   'ENSG00000163659',
		   'ENSG00000117298',
		   'ENSG00000166483',
		   'ENSG00000181026'};

highlight_labels = {'d', 'b', 'c', 'e', 'a', 'f'};

%id = '2013-08-30';
id = 'final';

profiles = load(['~/projects/pol2rnaseq/analyses/hmc_results/profiles/all_profiles_' id '.mat']);

%t_pred = (((0:100)/100*sqrt(1280)).^2)';
T_MIN = 300;
t_pred = profiles.t_pred - T_MIN;
realt = find(t_pred >= 0);

mygenes = importdata('../R/analysed_genes.txt');

[I, A, B] = intersect(mygenes, profiles.genes);

dataVals1 = profiles.mu(B, realt);
dataVals2 = profiles.mu(B, realt + length(profiles.t_pred));
dataVals1(dataVals1(:) < 0) = 0;

[~, Ahigh, Bhigh] = intersect(mygenes(A), highlight_genes);
assert(length(Bhigh) == length(highlight_genes))

I_pcomb = do_plot_heatmaps(dataVals1, dataVals2, t_pred(realt), ...
                           'plots/profiles', [], ...
			   struct('I', {Ahigh}, 'txt', {highlight_labels(Bhigh)}));

% -----------------------------------------------------------------
% DATA
% -----------------------------------------------------------------

datadir = '~/projects/pol2rnaseq/data/';

%load h3k4me3_series.mat
%load series_for_matti_ver3.mat
%load('pol2_for_matti_ver3.mat', 'bininfo');
%load('pol2_summaryseries_2012_09.mat');
load([datadir, 'bininfo_nov2014_corrected.mat'], 'bininfo');
load([datadir 'pol2_summaryseries_2014_11_19.mat']);
%r = load('rna_new_data4.mat');
%act = importdata('activeGenes_new.txt');
normfacts = importdata([datadir 'rna_norm_factors.txt']);
r = load([datadir 'info_gene_mean_var.mat']);

mygenes = importdata('../R/analysed_genes.txt');

[I, A, B] = intersect(ensg2int(r.geneID), bininfo(:, 5));

J = zeros(length(mygenes), 1);
for k=1:length(mygenes),
  JJ = strcmp(mygenes{k}, r.geneID(A));
  if any(JJ),
    J(k) = find(JJ);
  end
end

%interestinggenes = A; %find(r.pvals(A) < 0.1);
%interestinggenes = find(r.pvals(A) < 0.9);
%interestinggenes = find(sum(r.counts(A,:), 2) >= 1000);
interestinggenes_rna = A(sort(J(find(J))));
interestinggenes_pol2 = B(sort(J(find(J))));

dataVals1=pol2_summaryseries(interestinggenes_pol2,:);
dataVals2=r.mu(interestinggenes_rna,:) ./ repmat(normfacts', [length(interestinggenes_rna), 1]);

t_data = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280];
%t_plot = linspace(0, sqrt(1280), 101).^2;
t_plot = t_pred(realt);
[foo, J] = min(bsxfun(@(x, y) abs(x-y), t_data, t_plot(:)), [], 2);

dataVals1 = dataVals1(:, J);
dataVals2 = dataVals2(:, J);

do_plot_heatmaps(dataVals1, dataVals2, t_plot, 'plots/data', I_pcomb);
