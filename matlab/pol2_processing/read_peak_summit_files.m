function mappings=read_mappingfile(peakfilename,summitfilename);

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT'};

n_linestoskip = 0;

%------------------------------
% Initialize data cell array
%------------------------------
n_chromosomes=length(chromosomenames);
mappings=cell(n_chromosomes,7);

index_peakstart=1;
index_peakend=2;
index_peakscore=3;
index_summitstart=4;
index_summitend=5;
index_summitheight=6;
index_npeaks=7;
index_maxpeaks=8;
index_linenumber=9;


for i=1:n_chromosomes,
  peakstarts=int32([]);
  peakends=int32([]);
  peakscores=single([]);
  summitstarts=int32([]);
  summitends=int32([]);
  summitheights=single([]);
  linenumbers=int32([]);
  mappings{i,index_peakstart}=peakstarts;
  mappings{i,index_peakend}=peakends;
  mappings{i,index_peakscore}=peakscores;
  mappings{i,index_summitstart}=summitstarts;
  mappings{i,index_summitend}=summitends;
  mappings{i,index_summitheight}=summitheights;
  mappings{i,index_linenumber}=linenumbers;
  mappings{i,index_npeaks}=0;
  mappings{i,index_maxpeaks}=0;
end;


%------------------------------
% Check the files exist
%------------------------------
peaks_available = 1;
f1 = fopen(peakfilename,'r');
if f1==-1, 
  peaks_available=0; 
  fprintf(1,'Could not open peak file [%s] for reading\n',peakfilename);
end;

summits_available = 1;
f2 = fopen(summitfilename,'r');
if f2==-1, 
  summits_available=0; 
  fprintf(1,'Could not open summit file [%s] for reading\n',summitfilename);
end;


%------------------------------
% Skip first lines (headers)
%------------------------------
for i=1:n_linestoskip,
  tline1=fgets(f1);
end;
for i=1:n_linestoskip,
  tline2=fgets(f2);
end;


%------------------------------
% Read in the rest of the files
%------------------------------
chr_index=1;
peakstarts=mappings{chr_index,index_peakstart};
peakends=mappings{chr_index,index_peakend};
peakscores=mappings{chr_index,index_peakscore};
summitstarts=mappings{chr_index,index_summitstart};
summitends=mappings{chr_index,index_summitend};
summitheights=mappings{chr_index,index_summitheight};
linenumbers=mappings{chr_index,index_linenumber};
npeaks=mappings{chr_index,index_npeaks};
maxpeaks=mappings{chr_index,index_maxpeaks};

nlines=0;
tline1=fgets(f1);
tline2=fgets(f2);
lines_left=0;
if tline1 ~= -1, lines_left=1; end;
if tline2 ~= -1, lines_left=1; end;

while lines_left == 1,
  nlines=nlines+1;
  if mod(nlines,1000)==0, sprintf('Reading line %d\n',nlines), end;

  %------------------------------
  % parse the peak line, 6 columns of interest
  %------------------------------
  if tline1 ~= -1,

    linelength=length(tline1);
  
    % chromosome index for peak
    i1 = 1; 
    i2 = i1;
    while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 1
    line_chrindex1=-1;
    for k=1:n_chromosomes, 
      if strcmp(chromosomenames{k},tline1(i1:i2-1))==1,
        line_chrindex1=k;
        break;
      end;
    end;
    
    % peak start
    i1 = i2;
    while (double(tline1(i1)) == double(' ')) || (double(tline1(i1)) == 9), i1=i1+1; end; % start of column 2
    i2 = i1;
    while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 2
    line_peakstart = int32(str2double(tline1(i1:i2-1)));

    % peak end
    i1 = i2;
    while (double(tline1(i1)) == double(' ')) || (double(tline1(i1)) == 9), i1=i1+1; end; % start of column 3
    i2 = i1;
    while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 3
    line_peakend = int32(str2double(tline1(i1:i2-1)));

    % peak name, we do not store this
    i1 = i2;
    while (double(tline1(i1)) == double(' ')) || (double(tline1(i1)) == 9), i1=i1+1; end; % start of column 4
    i2 = i1;
    while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 4

    % peak score
    i1 = i2;
    while (i1 < linelength) && ((double(tline1(i1)) == double(' ')) || (double(tline1(i1)) == 9)), i1=i1+1; end; % start of column 5
    i2 = i1;
    while (i2 < linelength) && ((double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9)), i2=i2+1; end; % end of column 5
    line_peakscore = single(str2double(tline1(i1:i2-1)));

  end;

  
  if tline2 ~= -1,

    linelength=length(tline2);
  
    % chromosome index for summit
    i1 = 1; 
    i2 = i1;
    while (double(tline2(i2)) ~= double(' ')) && (double(tline2(i2)) ~= 9), i2=i2+1; end; % end of column 1
    line_chrindex2=-1;
    for k=1:n_chromosomes, 
      if strcmp(chromosomenames{k},tline2(i1:i2-1))==1,
        line_chrindex2=k;
        break;
      end;
    end;

    % summit start
    i1 = i2;
    while (double(tline2(i1)) == double(' ')) || (double(tline2(i1)) == 9), i1=i1+1; end; % start of column 2
    i2 = i1;
    while (double(tline2(i2)) ~= double(' ')) && (double(tline2(i2)) ~= 9), i2=i2+1; end; % end of column 2
    line_summitstart = int32(str2double(tline2(i1:i2-1)));

    % summit end
    i1 = i2;
    while (double(tline2(i1)) == double(' ')) || (double(tline2(i1)) == 9), i1=i1+1; end; % start of column 3
    i2 = i1;
    while (double(tline2(i2)) ~= double(' ')) && (double(tline2(i2)) ~= 9), i2=i2+1; end; % end of column 3
    line_summitend = int32(str2double(tline2(i1:i2-1)));

    % peak name, we do not store this
    i1 = i2;
    while (double(tline2(i1)) == double(' ')) || (double(tline2(i1)) == 9), i1=i1+1; end; % start of column 4
    i2 = i1;
    while (double(tline2(i2)) ~= double(' ')) && (double(tline2(i2)) ~= 9), i2=i2+1; end; % end of column 4

    % summit height
    i1 = i2;
    while (i1 < linelength) && ((double(tline2(i1)) == double(' ')) || (double(tline2(i1)) == 9)), i1=i1+1; end; % start of column 5
    i2 = i1;
    while (i2 < linelength) && ((double(tline2(i2)) ~= double(' ')) && (double(tline2(i2)) ~= 9)), i2=i2+1; end; % end of column 5
    line_summitheight = single(str2double(tline2(i1:i2-1)));
  
  end;
  
  % basic consistency check, chromosome indices should be the same for the peak and summit lines
  if (peaks_available == 1) && (summits_available == 1),
    if line_chrindex1 ~= line_chrindex2,
      fprintf(1,'Chromosome index inconsistency on line %d: %d (peak) vs. %d (summit)\n', nlines, chr_index1, chr_index2);
      return;
    end;
  end;

  if tline1 ~= -1,
    line_chrindex = line_chrindex1;
  end;
  if tline2 ~= -1,
    line_chrindex = line_chrindex2;
  end;
  
  % start working on the chromosome identified on the peak&summit lines, if we weren't already working on it
  if line_chrindex ~= chr_index,
    mappings{chr_index,index_peakstart}=peakstarts;
    mappings{chr_index,index_peakend}=peakends;
    mappings{chr_index,index_peakscore}=peakscores;
    mappings{chr_index,index_summitstart}=summitstarts;
    mappings{chr_index,index_summitend}=summitends;
    mappings{chr_index,index_summitheight}=summitheights;
    mappings{chr_index,index_linenumber}=linenumbers;
    mappings{chr_index,index_npeaks}=npeaks;
    mappings{chr_index,index_maxpeaks}=maxpeaks;
    
    chr_index = line_chrindex;
    peakstarts=mappings{chr_index,index_peakstart};
    peakends=mappings{chr_index,index_peakend};
    peakscores=mappings{chr_index,index_peakscore};
    summitstarts=mappings{chr_index,index_summitstart};
    summitends=mappings{chr_index,index_summitend};
    summitheights=mappings{chr_index,index_summitheight};
    linenumbers=mappings{chr_index,index_linenumber};
    npeaks=mappings{chr_index,index_npeaks};
    maxpeaks=mappings{chr_index,index_maxpeaks};    
  end;

  % add the mapped read to the data of the chromosome
  npeaks=npeaks+1;    
  if (npeaks>maxpeaks), 
    maxpeaks=npeaks*2;
    peakstarts=[peakstarts zeros(1,maxpeaks-npeaks+1,'int32')];
    peakends=[peakends zeros(1,maxpeaks-npeaks+1,'int32')];
    peakscores=[peakscores zeros(1,maxpeaks-npeaks+1,'single')];
    summitstarts=[summitstarts zeros(1,maxpeaks-npeaks+1,'int32')];
    summitends=[summitends zeros(1,maxpeaks-npeaks+1,'int32')];
    summitheights=[summitheights zeros(1,maxpeaks-npeaks+1,'single')];
    linenumbers=[linenumbers zeros(1,maxpeaks-npeaks+1,'int32')];
  end;
  peakstarts(npeaks)=line_peakstart;
  peakends(npeaks)=line_peakend;
  peakscores(npeaks)=line_peakscore;
  summitstarts(npeaks)=line_summitstart;
  summitends(npeaks)=line_summitend;
  summitheights(npeaks)=line_summitheight;
  linenumbers(npeaks)=nlines;  % line number where this read occurred in the file

  % read in the next lines
  tline1=fgets(f1);
  tline2=fgets(f2);
  lines_left=0;
  if tline1 ~= -1, lines_left=1; end;
  if tline2 ~= -1, lines_left=1; end;

  % debug early stopping point
  % if nlines > 1000, tline1=-1,end;
end;
  

% store the data of the last chromosome we were working on
mappings{chr_index,index_peakstart}=peakstarts;
mappings{chr_index,index_peakend}=peakends;
mappings{chr_index,index_peakscore}=peakscores;
mappings{chr_index,index_summitstart}=summitstarts;
mappings{chr_index,index_summitend}=summitends;
mappings{chr_index,index_summitheight}=summitheights;
mappings{chr_index,index_linenumber}=linenumbers;
mappings{chr_index,index_npeaks}=npeaks;
mappings{chr_index,index_maxpeaks}=maxpeaks;

% for all chromosomes, discard empty data
for chr_index=1:n_chromosomes,
  npeaks=mappings{chr_index,index_npeaks};
  mappings{chr_index,index_maxpeaks}=npeaks;
  mappings{chr_index,index_peakstart}=mappings{chr_index,index_peakstart}(1:npeaks);
  mappings{chr_index,index_peakend}=mappings{chr_index,index_peakend}(1:npeaks);
  mappings{chr_index,index_peakscore}=mappings{chr_index,index_peakscore}(1:npeaks);
  mappings{chr_index,index_summitstart}=mappings{chr_index,index_summitstart}(1:npeaks);
  mappings{chr_index,index_summitend}=mappings{chr_index,index_summitend}(1:npeaks);
  mappings{chr_index,index_summitheight}=mappings{chr_index,index_summitheight}(1:npeaks);  
  mappings{chr_index,index_linenumber}=mappings{chr_index,index_linenumber}(1:npeaks);
end;

fclose(f1);
fclose(f2);



