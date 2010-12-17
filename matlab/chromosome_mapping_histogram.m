function [bins,binsstart]=chromosome_mapping_histogram(mappings,chr_index,interval,fragmentlength);

readstarts_index=1;
readends_index=2;
readstrands_index=3;
readscores_index=4;
nreads_index=5;
maxreads_index=6;
linenumbers_index=7;

binsstart=floor(double((min(mappings{chr_index,readstarts_index}))-fragmentlength)/interval)*interval;
binsend=floor(double((max(mappings{chr_index,readends_index}))+fragmentlength)/interval)*interval;

nbins = (binsend-binsstart)/interval +1;

bins = zeros(nbins,1);

nreads=length(mappings{chr_index,readstarts_index});
for i=1:nreads,
  if mod(i,10000)==0, i, end;
  readstart=mappings{chr_index,readstarts_index}(i);
  readend=mappings{chr_index,readends_index}(i);
  readscore=mappings{chr_index,readscores_index}(i);
  readstrand=mappings{chr_index,readstrands_index}(i);
  
  if (readstrand>0),
    % shift
    readstart=readstart+fragmentlength/2;
    readend=readend+fragmentlength/2;
  
    % extend
    readstart=readstart-fragmentlength/2;
    readend=readstart+fragmentlength;  
  end;
  if (readstrand<0),
    % shift
    readstart=readstart-fragmentlength/2;
    readend=readend-fragmentlength/2;

    % extend
    readend=readend+fragmentlength/2;  
    readstart=readend-fragmentlength;
  end;
  
  firstbin=floor(double((readstart-binsstart))/interval);
  lastbin=floor(double((readend-binsstart))/interval);

  
  

  if lastbin > firstbin,  
    firstprop=binsstart+interval*firstbin + interval - readstart;
    lastprop=readend - binsstart-interval*lastbin + 1;

    %[binsstart binsend readstart readend firstbin firstprop lastbin lastprop]

    if firstprop < 0 || lastprop < 0,
      [binsstart binsend readstart readend firstbin firstprop lastbin lastprop]
      pause
    end;
    
    %[double(firstprop) double(readscore)]
    bins(firstbin+1)=bins(firstbin+1)+firstprop*double(readscore);
    bins(lastbin+1)=bins(lastbin+1)+lastprop*double(readscore);
    bins(firstbin+2:lastbin)=bins(firstbin+2:lastbin)+interval*double(readscore);
  else
    bins(firstbin+1)=bins(firstbin+1)+(readend-readstart+1)*double(readscore);
  end;
end;





