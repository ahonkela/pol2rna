function [mappings,duplicates_removed] = discard_duplicates(mappings, max_duplicates, already_sorted);

readstarts_index=1;
readends_index=2;
readstrands_index=3;
readscores_index=4;
nreads_index=5;
maxreads_index=6;
linenumbers_index=7;

n_chromosomes=size(mappings,1);
duplicates_removed=zeros(n_chromosomes,1);

for chr_index=1:size(mappings,1),
  chr_index
  if already_sorted==1,
    I1=[1:mappings{chr_index,nreads_index}];
  else    
    [y,I1] = sortrows([mappings{chr_index,readstarts_index} mappings{chr_index,readends_index}],[1 2]);
  end;
  I2=zeros(length(I1),1);
  I2(1)=I1(1);
  
  duplicates=1;
  readstart=mappings{chr_index,readstarts_index}(I1(1));
  readend=mappings{chr_index,readstarts_index}(I1(1));
  nreads_final=1;
  for j=2:length(I1),
    if mod(j,500000)==0,
      j
    end;
    
    tempreadstart=mappings{chr_index,readstarts_index}(I1(j));
    tempreadend=mappings{chr_index,readstarts_index}(I1(j));

    if (tempreadstart==readstart) && (tempreadend==readend),
      duplicates=duplicates+1;
    else
      duplicates=1;
      readstart=tempreadstart;
      readend=tempreadend;
    end;
    if duplicates <= max_duplicates,
      nreads_final=nreads_final+1;
      I2(nreads_final)=I1(j);
    else
      duplicates_removed(chr_index)=duplicates_removed(chr_index)+1;
    end;
  end;

  I2=I2(1:nreads_final);
  
  mappings{chr_index,readstarts_index}=mappings{chr_index,readstarts_index}(I2);
  mappings{chr_index,readends_index}=mappings{chr_index,readends_index}(I2);
  mappings{chr_index,readstrands_index}=mappings{chr_index,readstrands_index}(I2);
  mappings{chr_index,readscores_index}=mappings{chr_index,readscores_index}(I2);
  mappings{chr_index,linenumbers_index}=mappings{chr_index,linenumbers_index}(I2);
end;

