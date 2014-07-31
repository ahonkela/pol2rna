expdata = load('~/projects/pol2rnaseq/data/sorted_tr_rpkm.mat');

N_tr = sum(cellfun('length', expdata.t0000_sorted_tr_names));
tr_exp = zeros(N_tr, 10);
tr_names = cell(N_tr, 1);

times = [0, 5, 10, 20, 40, 80, 160, 320, 640, 1280];
namefields = cell(length(times), 1);
rpkmfields = cell(length(times), 1);
for k=1:length(times),
  namefields{k} = sprintf('t%04d_sorted_tr_names', times(k));
  rpkmfields{k} = sprintf('t%04d_sorted_tr_rpkm', times(k));
end

begind = 0;
for k=1:length(expdata.genes),
  if rem(k, 1000)==0,
    fprintf('Doing gene %d/%d\n', k, length(expdata.genes));
  end
  N = length(expdata.(namefields{1}){k});
  trs = sort(expdata.(namefields{1}){k});
  for m=1:N,
    tr_names{begind + m} = sprintf('%s.%s', expdata.genes{k}, trs{m});
  end
  for l=1:length(times),
    [~, I] = sort(expdata.(namefields{l}){k});
    tr_exp(begind+(1:N),l) = expdata.(rpkmfields{l}){k}(I);
  end
  begind = begind + N;
end

fp = fopen('~/projects/pol2rnaseq/data/transcript_expression.txt', 'w');
for k=1:N_tr,
  fprintf(fp, '%s %f %f %f %f %f %f %f %f %f %f\n', tr_names{k}, tr_exp(k,:));
end
fclose(fp);
