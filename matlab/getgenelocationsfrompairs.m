% Extract the locations of genes involved in pol2-rna pairs
% Jaakko Peltonen, Jan 12-, 2011.

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT'};
n_chromosomes=length(chromosomenames);

geneinfo=cell(length(pairs),1);
for i=1:length(pairs),
  genename=pairs{i}{1,3};
  genechromosome=pairs{i}{1,4};
  genestart=pairs{i}{1,5};
  geneend=pairs{i}{1,6};
  geneinfotemp=cell(5,1);
  geneinfotemp{1}=genename;
  geneinfotemp{2}=genechromosome;
  line_chrindex=-1;
  for k=1:n_chromosomes,
    if strcmp(genechromosome,chromosomenames{k})==1,
      line_chrindex=k;
      break;
    end;
  end;
  geneinfotemp{3}=line_chrindex;
  geneinfotemp{4}=genestart;
  geneinfotemp{5}=geneend;

  geneinfo{i}=geneinfotemp;
end;

uniques = ones(size(pairs,2),1);
nuniques=0;

for i=1:size(pairs,2),
  if mod(i,100)==0, i, end;

duplicate_found=0;
  for j=1:nuniques,
    if (geneinfo{i}{4}==geneinfo{uniques(j)}{4}) && (geneinfo{i}{5}==geneinfo{uniques(j)}{5}),
      duplicate_found=1;
      break;
    end;
  end;
  if duplicate_found==0,
    nuniques=nuniques+1;
    uniques(nuniques)=i;
  end;    
end;

geneinfo={geneinfo{uniques(1:nuniques)}};

