function write_mappingfile(filename,readdata);

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT','rRNA'};
n_chromosomes=length(chromosomenames);

f = fopen(filename,'w');

for k=1:size(readdata,1),
  n_genes=length(readdata{k,1});
  for l=1:n_genes,
    chrname=chromosomenames{k};
    readstart=readdata{k,1}(l);
    readend=readdata{k,2}(l);
    if readdata{k,3}(l)==1,
      readstrand='+';
    else
      readstrand='-';
    end;
    readscore=readdata{k,4}(l);
    fprintf(f,'%s\t%d\t%d\tn\t%f\t%s\n',chrname,readstart,readend,readscore,readstrand);
  end;
end;

fclose(f);
