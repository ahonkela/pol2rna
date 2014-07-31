function bins=computepol2activityovergenes(mappings,chr_index,binstarts,binends,fragmentlength,n_allowed_duplicates);

readstarts_index=1;
readends_index=2;
readstrands_index=3;
readscores_index=4;
nreads_index=5;
maxreads_index=6;
linenumbers_index=7;

nbins=length(binstarts);

bins = zeros(nbins,1);

nreads=length(mappings{chr_index,readstarts_index});
readstart=-1;
readend=-1;

n_duplicatesfound=0;

for i=1:nreads,
  if mod(i,50)==0, i, end;
  
  tempreadstart=mappings{chr_index,readstarts_index}(i);
  tempreadend=mappings{chr_index,readends_index}(i);
  
  %-------------------------------------
  % check for a duplicate (assume reads are sorted so duplicates
  % are listed one after another)
  %-------------------------------------
  if (tempreadstart==readstart) && (tempreadend==readend),
    n_duplicatesfound=n_duplicatesfound+1;
  else
    n_duplicatesfound=1;
    readstart=tempreadstart;
    readend=tempreadend;
  end;
  
  if n_duplicatesfound <= n_allowed_duplicates,
  
    readscore=mappings{chr_index,readscores_index}(i);
    readstrand=mappings{chr_index,readstrands_index}(i);

    %-------------------------------------
    % shift and extend the read using the desired fragment length
    %-------------------------------------
    if fragmentlength > 0,  
      if (readstrand>0),
	% shift
	modified_readstart=readstart+fragmentlength/2;
	modified_readend=readend+fragmentlength/2;
  
	% extend
	modified_readstart=modified_readstart-fragmentlength/2;
	modified_readend=modified_readstart+fragmentlength;  
      end;
      if (readstrand<0),
	% shift
	modified_readstart=readstart-fragmentlength/2;
	modified_readend=readend-fragmentlength/2;
	
	% extend
	modified_readend=modified_readend+fragmentlength/2;  
	modified_readstart=modified_readend-fragmentlength;
      end;
    end;

    %-------------------------------------
    % Compute assignment to bins.
    % Assume the bins have been sorted according to start index.
    % Assume that bins may overlap. A read that falls into multiple bins
    % is currently counted separately for each bin (i.e. the
    % influence of the read is multiplied).
    %-------------------------------------

    % find latest bin (greatest index) whose start is not after the read
    k0=nbins; while (k0>=1) && (binstarts(k0)>modified_readend),k0=k0-1;end;
    
    for k=k0:-1:1,
      % Since we know the end of the read is at or after the start
      % of the bin, all we need to check is that the start of the 
      % read has not happened after the bin (i.e. check that the
      % whole read is not to the right of the bin).
      if modified_readstart>binends(k),
	continue;
      end;

      
      
      %-------------------------------------
      % Compute overlapping portion between the read and the bin.
      %-------------------------------------
      if modified_readend>=binends(k), 
	tempend=binends(k); 
      else 
	tempend=modified_readend;
      end;
      if modified_readstart<=binstarts(k),
	tempstart=binstarts(k);
      else
	tempstart=modified_readstart;
      end;

      if tempstart>tempend,
	fprintf(1,'problem: bin %d, %d, read %d, %d, overlap %d\n', ...
		   binstarts(k),binends(k),modified_readstart, ...
		modified_readend,tempend-tempstart+1);
	pause;
      end;      
      
      %-------------------------------------
      % Assign the read to the bin.
      % Currently, amount added to the bin depends on length of
      % overlap and on score of the read.      
      %-------------------------------------
      bins(k)=bins(k)+(tempend-tempstart+1)*double(readscore);
    end;  
    
  end;  
end;

