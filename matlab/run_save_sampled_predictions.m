pol2rnaseqToolboxes;

fid = fopen('~/github/pol2rna/matlab/file_predictions_config.txt');
conf = textscan(fid, '%s');
fclose(fid);
conf = conf{1};

g = importdata(['~/github/pol2rna/matlab/', conf{1}]);
myI = mybase:mymod:length(g);

FILESTEM = conf{2};
SAMPLEDIR = conf{3};  
SQRTTIME = str2num(conf{4});
FILESPEC = conf{5};

for k=myI,
  fprintf('Running gene %d/%d: %s\n', find(k==myI), length(myI), g{k});
  save_sampled_predictions_file(g{k}, SAMPLEDIR, FILESTEM, FILESPEC, SQRTTIME);
end
