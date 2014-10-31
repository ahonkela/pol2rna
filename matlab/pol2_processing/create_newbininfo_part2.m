addpath /share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/pol2rnaseq/matlab

binfile = read_stringfile('/share/mi/workspace/jtpelto/synergy/synergy_data/analyses/new_bininfo/genes_dec2012_all.txt',[' ' 10 13]);

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};

%   bininfo(i,1)=line_chrindex; % chromosome index
%   bininfo(i,2)=str2double(genetemp{3}(2:end-1)); % start location
%   bininfo(i,3)=str2double(genetemp{4}(2:end-1)); % end location
%   bininfo(i,4)=i; % line number in an original file, never mind this
%   bininfo(i,5)=geneids(i);  % ENSEMBL id of the gene
%   bininfo(i,6)=genestrands(i); % strand identifier (+1 or -1)

bininfo = nan*ones(length(binfile),6);
for k=1:length(binfile),
  if mod(k,1000)==0,
    k
  end;
  % line number in original file, ignore this
  bininfo(k,4)=k;   
  % ENSEMBL id
  bininfo(k,5)=str2double(binfile{k}{2}(6:end-1));  
  % chromosome id
  chrtemp=binfile{k}{3}(2:end-1);
  bininfo(k,1)=-1;
  for l=1:length(chromosomenames),
    if strcmp(chrtemp,chromosomenames{l})==1,
      bininfo(k,1)=l;
    end;
  end;
  % strand id
  bininfo(k,6)=0;
  strandtemp=binfile{k}{4}(2:end-1);
  if strandtemp(1)=='+',
    bininfo(k,6)=1;
  end;
  if strandtemp(1)=='-',
    bininfo(k,6)=-1;
  end;
  % start position
  bininfo(k,2)=str2double(binfile{k}{5}(2:end-1));
  % end position
  bininfo(k,3)=str2double(binfile{k}{6}(2:end-1));
end;

[ytemp,Itemp]=sort(bininfo(:,2));
bininfo_sorted=bininfo(Itemp,:);
bininfo_unsorted = bininfo;

bininfo = bininfo_sorted;

save bininfo_dec2012.mat bininfo bininfo_unsorted bininfo_sorted -mat
% save bininfo_dec2012.txt bininfo -ASCII
f=fopen('bininfo_dec2012.txt','w');
for k=1:size(bininfo,1),
  fprintf(f,'%d %d %d %d %d %d\n',bininfo(k,1),bininfo(k,2),bininfo(k,3),bininfo(k,4),bininfo(k,5),bininfo(k,6));
end;
fclose(f);


