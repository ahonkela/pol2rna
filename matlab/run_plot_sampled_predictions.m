pol2rnaseqToolboxes;

if ~exist('config', 'var'),
  config = 'file_predictions_config.txt';
end

fid = fopen(['~/github/pol2rna/matlab/' config]);
conf = textscan(fid, '%s');
fclose(fid);
conf = conf{1};

g = importdata(['~/github/pol2rna/matlab/', conf{1}]);
myI = mybase:mymod:length(g);

if length(conf) > 3,
  SQRTTIME = str2num(conf{4});
else
  SQRTTIME = 1;
end

FILESTEM = conf{2};
SAMPLEDIR = conf{3};  
FILESPEC = conf{5};

if length(conf) > 5,
  PLOTPREMRNA = str2num(conf{6});
else
  PLOTPREMRNA = 0;
end

if length(conf) > 6,
  PREMRNAMODEL = str2num(conf{7});
else
  PREMRNAMODEL = 0;
end

for k=myI,
  fprintf('Running gene %d/%d: %s\n', find(k==myI), length(myI), g{k});
  plot_sampled_predictions_file(g{k}, SAMPLEDIR, FILESTEM, 'png', SQRTTIME, FILESPEC, PLOTPREMRNA, PREMRNAMODEL);
end
