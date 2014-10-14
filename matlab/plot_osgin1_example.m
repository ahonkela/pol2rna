FONTSIZE=6;

datadir = '~/projects/pol2rnaseq/data/';
normfacts = importdata([datadir, 'rna_norm_factors.txt']);
r = load([datadir, 'info_gene_mean_var.mat']);

counts = importdata('../python/results/OSGIN1_data.txt');

rpm = (counts(:, 1:2) ./ repmat((counts(:, 3) .* normfacts)/1e6, [1, 2]))';
rpm(1,:) = rpm(1,:) - rpm(2,:);

OSGIN1_ind = find(strcmp(r.gene_name, 'OSGIN1'));

bitseq = r.mu(OSGIN1_ind, :) ./ normfacts';

t = [0, 5, 10, 20, 40, 80, linspace(160, 320, 4)];
t_tick = [0, 40, 80, linspace(160, 320, 4)];
t_label = [0, 40, 80, 160, 320, 640, 1280];

r2 = load('/home/ahonkela/projects/pol2rnaseq/data/sorted_tr_rpkm.mat');
keys = {'t0000', 't0005', 't0010', 't0020', 't0040', 't0080', 't0160', 't0320', 't0640', 't1280'};
short_share = zeros(10, 1);
for k=1:10,
  shortI = find(strcmp(r2.([keys{k}, '_sorted_tr_names']){OSGIN1_ind}, 'ENST00000563543'));
  my_rpkms = r2.([keys{k}, '_sorted_tr_rpkm']){OSGIN1_ind};
  short_share(k) = my_rpkms(shortI) ./ sum(my_rpkms);
end

subplot(1, 2, 1)
plot(t, bitseq ./ sqrt(mean(bitseq.^2)), 'bo-');
set(gca, 'FontSize', FONTSIZE)
hold on
plot(t, rpm(1,:) ./ sqrt(mean(rpm(1,:).^2)), 'ro-');
%plot(t, rpm(2,:) ./ sqrt(sum(rpm(2,:).^2)), 'go-');
hold off
axis([0 320 0 2])
legend('Transcript-based', 'Count-based', 'location', 'SouthEast')
set(gca, 'XTick', t_tick);
set(gca, 'XTickLabel', t_label)
xlabel('t (min)')
ylabel('Expression (rescaled)')

subplot(1, 2, 2)
plot(t, short_share, 'kx-');
set(gca, 'FontSize', FONTSIZE)
axis([0 320 0 1])
set(gca, 'XTick', t_tick);
set(gca, 'XTickLabel', t_label)
xlabel('t (min)')
ylabel('ENST00000563543 relative expression')

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0, 0, 11.7, 6.0])
print('-depsc2', 'plots/osgin1_example.eps');

