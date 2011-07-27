function mappings=read_mappingfile(filename);

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT'};

n_linestoskip = 4;

%------------------------------
% Initialize data cell array
%------------------------------
n_chromosomes=length(chromosomenames);
mappings=cell(n_chromosomes,7);

index_readstart=1;
index_readend=2;
index_readstrand=3;
index_readscore=4;
index_nreads=5;
index_maxreads=6;
index_linenumber=7;

for i=1:n_chromosomes,
  readstarts=int32([]);
  readends=int32([]);
  readstrands=int8([]);
  readscores=single([]);
  linenumbers=int32([]);
  mappings{i,index_readstart}=readstarts;
  mappings{i,index_readend}=readends;
  mappings{i,index_readstrand}=readstrands;
  mappings{i,index_readscore}=readscores;
  mappings{i,index_linenumber}=linenumbers;
  mappings{i,index_nreads}=0;
  mappings{i,index_maxreads}=0;
end;


%------------------------------
% Check the file exists
%------------------------------
f1 = fopen(filename,'r');
if f1==-1, 
  mappings_available=0; 
  fprintf(1,'Could not open file [%s] for reading\n',filename);
  return;
end;


%------------------------------
% Skip first lines (header)
%------------------------------
for i=1:n_linestoskip,
  tline1=fgets(f1);
end;


%------------------------------
% Read in the rest of the file
%------------------------------
chr_index=1;
readstarts=mappings{chr_index,index_readstart};
readends=mappings{chr_index,index_readend};
readstrands=mappings{chr_index,index_readstrand};
readscores=mappings{chr_index,index_readscore};
linenumbers=mappings{chr_index,index_linenumber};
nreads=mappings{chr_index,index_nreads};
maxreads=mappings{chr_index,index_maxreads};

nlines=0;
tline1=fgets(f1);
while tline1 ~= -1,
  nlines=nlines+1;
  if mod(nlines,1000)==0, sprintf('Reading line %d\n',nlines), end;

  %------------------------------
  % parse the line, 6 columns of interest
  %------------------------------

  % chromosome index
  i1 = 1; 
  i2 = i1;
  while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 1
  line_chrindex=-1;
  for k=1:n_chromosomes, 
    if strcmp(chromosomenames{k},tline1(i1:i2-1))==1,
      line_chrindex=k;
      break;
    end;
  end;

  % mapping start
  i1 = i2;
  while (double(tline1(i1)) == double(' ')) || (double(tline1(i1)) == 9), i1=i1+1; end; % start of column 2
  i2 = i1;
  while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 2
  line_readstart = int32(str2double(tline1(i1:i2-1)));

  % mapping end
  i1 = i2;
  while (double(tline1(i1)) == double(' ')) || (double(tline1(i1)) == 9), i1=i1+1; end; % start of column 3
  i2 = i1;
  while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 3
  line_readend = int32(str2double(tline1(i1:i2-1)));

  % read name, we do not store this
  i1 = i2;
  while (double(tline1(i1)) == double(' ')) || (double(tline1(i1)) == 9), i1=i1+1; end; % start of column 4
  i2 = i1;
  while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 4

  % mapping score
  i1 = i2;
  while (double(tline1(i1)) == double(' ')) || (double(tline1(i1)) == 9), i1=i1+1; end; % start of column 5
  i2 = i1;
  while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 5
  line_readscore = single(str2double(tline1(i1:i2-1)));

  % strand identifier, '+' or '-'
  i1 = i2;
  while (double(tline1(i1)) == double(' ')) || (double(tline1(i1)) == 9), i1=i1+1; end; % start of column 2
  i2 = i1;
  while (double(tline1(i2)) ~= double(' ')) && (double(tline1(i2)) ~= 9), i2=i2+1; end; % end of column 2
  if tline1(i1) == '+', 
    line_readstrand = int8(1); 
  elseif tline1(i1) == '-', 
    line_readstrand = int8(-1); 
  else 
    line_readstrand = int8(0); 
  end;

  % start working on the chromosome identified on the line, if we weren't already working on it
  if line_chrindex ~= chr_index,
    mappings{chr_index,index_readstart}=readstarts;
    mappings{chr_index,index_readend}=readends;
    mappings{chr_index,index_readstrand}=readstrands;
    mappings{chr_index,index_readscore}=readscores;
    mappings{chr_index,index_linenumber}=linenumbers;
    mappings{chr_index,index_nreads}=nreads;
    mappings{chr_index,index_maxreads}=maxreads;
    
    chr_index = line_chrindex;
    readstarts=mappings{chr_index,index_readstart};
    readends=mappings{chr_index,index_readend};
    readstrands=mappings{chr_index,index_readstrand};
    readscores=mappings{chr_index,index_readscore};
    linenumbers=mappings{chr_index,index_linenumber};
    nreads=mappings{chr_index,index_nreads};
    maxreads=mappings{chr_index,index_maxreads};
  end;

  % add the mapped read to the data of the chromosome
  nreads=nreads+1;    
  if (nreads>maxreads), 
    maxreads=nreads*2;
    readstarts=[readstarts zeros(1,maxreads-nreads+1,'int32')];
    readends=[readends zeros(1,maxreads-nreads+1,'int32')];
    readstrands=[readstrands zeros(1,maxreads-nreads+1,'int8')];
    readscores=[readscores zeros(1,maxreads-nreads+1,'single')];
    linenumbers=[linenumbers zeros(1,maxreads-nreads+1,'int32')];
  end;
  readstarts(nreads)=line_readstart;
  readends(nreads)=line_readend;
  readstrands(nreads)=line_readstrand;
  readscores(nreads)=line_readscore;
  linenumbers(nreads)=nlines;  % line number where this read occurred in the file

  % read in the next line
  tline1=fgets(f1);

  % debug early stopping point
  % if nlines > 1000, tline1=-1,end;
end;
  

% store the data of the last chromosome we were working on
mappings{chr_index,index_readstart}=readstarts;
mappings{chr_index,index_readend}=readends;
mappings{chr_index,index_readstrand}=readstrands;
mappings{chr_index,index_readscore}=readscores;
mappings{chr_index,index_linenumber}=linenumbers;
mappings{chr_index,index_nreads}=nreads;
mappings{chr_index,index_maxreads}=maxreads;

% for all chromosomes, discard empty data
for chr_index=1:n_chromosomes,
  nreads=mappings{chr_index,index_nreads};
  mappings{chr_index,index_maxreads}=nreads;
  mappings{chr_index,index_readstart}=mappings{chr_index,index_readstart}(1:nreads);
  mappings{chr_index,index_readend}=mappings{chr_index,index_readend}(1:nreads);
  mappings{chr_index,index_readstrand}=mappings{chr_index,index_readstrand}(1:nreads);
  mappings{chr_index,index_readscore}=mappings{chr_index,index_readscore}(1:nreads);
  mappings{chr_index,index_linenumber}=mappings{chr_index,index_linenumber}(1:nreads);
end;

fclose(f1);



