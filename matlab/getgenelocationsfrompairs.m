% Extract the locations of genes involved in pol2-rna pairs
% Jaakko Peltonen, Jan 12 2011.

geneinfo=cell(length(pairs),4);
for i=1:length(pairs),
  genename=pairs{i}{1,3};
  genechromosome=pairs{i}{1,4};
  genestart=pairs{i}{1,5};
  geneend=pairs{i}{1,6};
  geneinfo{i,1}=genename;
  geneinfo{i,2}=genechromosome;
  geneinfo{i,3}=genestart;
  geneinfo{i,4}=geneend;
end;

uniques = ones(size(pairs,2),1);
nuniques=0;

for i=1:size(pairs,2),
  if mod(i,100)==0, i, end;

duplicate_found=0;
  for j=1:nuniques,
    if (geneinfo{i,3}==geneinfo{uniques(j),3}) && (geneinfo{i,4}==geneinfo{uniques(j),4}),
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

