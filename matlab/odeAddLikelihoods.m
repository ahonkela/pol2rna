function odeAddLikelihoods(fnames),

pol2rnaseqToolboxes;

d = load('ode_data_2015-05-07.mat');

for k=1:length(fnames),
  fname = fnames{k};
  [path, name, ext] = fileparts(fname);
  t = strsplit(name, '_');
  gene = t{1};

  r = load(fname);
  if isfield(r, 'll'),
    fprintf('file %s already contains likelihoods, skipping...\n', fname);
    continue;
  end
  fprintf('running file %s\n', fname);
  
  myindex = find(strcmp(gene, d.gene_name));
  ll = zeros(size(r.samples, 1), 1);
  for l = 1:length(ll),
    if mod(l, 10) == 0,
      fprintf('sample %d/%d\n', l, length(ll));
    end
    ll(l) = odeLikelihood(d.dataVals1(:, myindex), d.dataVals2(:, myindex), ...
                          d.rnaVars(:, myindex), d.timevector, ...
                          r.samples(l, :), 0, 0);
  end
  gene_name = r.gene_name;
  gene_index = r.gene_index;
  samples = r.samples;
  accepts = r.accepts;
  safeSave(fname, gene_name, gene_index, samples, accepts, ll);
end
