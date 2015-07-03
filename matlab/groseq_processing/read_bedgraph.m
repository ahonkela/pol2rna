function beddata=read_bedgraph(filename);

%
% Read a BEDGRAPH file into a numeric format
% when the original format has lines of this form:
% chr1    13487   13537   -1
%

f=fopen(filename,'r');

% This maximum file length is actually not a maximum, it simply
% helps Matlab to allocate a big chunk at the start. If lines
% exceed the maximum the code slows down due to reallocations of
% the beddata array.

maxlines=3260000;
beddata=zeros(maxlines,4);


% Read in each line, and compare the chromosome string to known chromosomes

iline=0;
tline=fgets(f);
knownchromosomes={'chr1','chr2','chr3','chr4','chr5','chr6','chr7', ...
                  'chr8','chr9','chr10','chr11','chr12','chr13', ...
                  'chr14','chr15','chr16','chr17','chr18','chr19', ...
                  'chr20','chr21','chr22','chrX','chrY'};
nchromosomes=length(knownchromosomes);
%maxlines=10;
while tline~=-1,
    iline=iline+1; 
    if mod(iline,10000)==0,
        iline
    end;
    
    % The last three numbers are the numeric values, the ones
    % before that are part of the chromosome name
    linenumbers=sscanf(tline,'%s %d %d %d');
    chromosome=char(linenumbers(1:end-3))';
    for k=1:nchromosomes,
        if strcmp(chromosome,knownchromosomes{k})==1,
            beddata(iline,1)=k;
            k=nchromosomes+1;
        end;        
    end;       
    beddata(iline,2:end)=linenumbers(end-2:end);
    tline=fgets(f);
    %if iline>maxlines, break; end;
end;
fclose(f);

% Remove unneeded lines from the array
%
beddata=beddata(1:iline,:);

