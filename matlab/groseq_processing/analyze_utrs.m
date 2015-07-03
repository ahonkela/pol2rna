% Step 1:
% Read in gene and UTR3 information generated from Ensembl Biomart
%--------------------------------

headers=read_stringfile('ensemblv68_utrs_headers_v2.txt',[10 13],['|']);
headerlengths=zeros(length(headers),1);
for k=1:length(headers),
    headerlengths(k)=length(headers{k});
end;


% Step 2:
% Parse text-format gene and UTR3 information into numeric information
%--------------------------------

ngenes=length(headers);
utrlocations=zeros(ngenes,7);
knownchromosomes={'1','2','3','4','5','6','7', ...
                  '8','9','10','11','12','13', ...
                  '14','15','16','17','18','19', ...
                  '20','21','22','X','Y'};
nchromosomes=length(knownchromosomes);

for k=1:ngenes,
    % ensembl id
    if strncmp(headers{k}{1},'>ENSG',5)==1,
        utrlocations(k,1)=str2double(headers{k}{1}(6:end));
    else
        utrlocations(k,1)=nan;
    end;    
    % utr start - may contain several colon-separated values, take
    % the minimum
    if length(headers{k}{2})>0,
        templocations=sscanf(headers{k}{2},'%d;');
        utrlocations(k,2)=min(templocations);
    else
        utrlocations(k,2)=nan;
    end;
    % utr end - may contain several colon-separated values, take
    % the maximum
    if length(headers{k}{3})>0,
        templocations=sscanf(headers{k}{3},'%d;');
        utrlocations(k,3)=max(templocations);
    else
        utrlocations(k,3)=nan;
    end;
    % gene start
    if length(headers{k}{4})>0,
        utrlocations(k,5)=str2double(headers{k}{4});
    else
        utrlocations(k,5)=nan;
    end;
    % gene end
    if length(headers{k}{5})>0,
        utrlocations(k,6)=str2double(headers{k}{5});
    else
        utrlocations(k,6)=nan;
    end;    
    % strand (1 or -1)
    if length(headers{k}{6})>0,
        utrlocations(k,4)=str2double(headers{k}{6});
    else
        utrlocations(k,4)=nan;
    end;
    % chromosome
    utrlocations(k,7)=nan;
    chromosomename=headers{k}{7};
    for m=1:nchromosomes,
        if strcmp(chromosomename,knownchromosomes{m})==1,
            utrlocations(k,7)=m;
            m=nchromosomes+1;
        end;
    end;   
end;


% Step 3:
% There may be multiple lines from different transcripts for each
% gene. Combine them to get the maximum gene area and maximum UTR3 area.
%--------------------------------

% Find genes that have valid ensembl IDs
I=find(isnan(utrlocations(:,1))==0);        
% Find the unique ensembl IDs (unique genes among the many transcripts)
uniquegeneids=unique(utrlocations(I,1));
uniquegenes=zeros(length(uniquegeneids),7);
% Find the maximum gene area and maximum UTR3 area for each gene
for k=1:length(uniquegeneids),
    I=find(utrlocations(:,1)==uniquegeneids(k));
    uniquegenes(k,1)=uniquegeneids(k);
    uniquegenes(k,4)=utrlocations(I(1),4); % strand
    uniquegenes(k,7)=utrlocations(I(1),7); % chromosome
    uniquegenes(k,2)=inf; % utr3 start
    uniquegenes(k,3)=-inf; % utr3 end
    uniquegenes(k,5)=inf; % gene start
    uniquegenes(k,6)=-inf; % gene end
    for m=1:length(I),
        if utrlocations(I(m),2) < uniquegenes(k,2),
            uniquegenes(k,2)= utrlocations(I(m),2);
        end;        
        if utrlocations(I(m),3) > uniquegenes(k,3),
            uniquegenes(k,3)= utrlocations(I(m),3);
        end;        
        if utrlocations(I(m),5) < uniquegenes(k,5),
            uniquegenes(k,5)= utrlocations(I(m),5);
        end;        
        if utrlocations(I(m),6) > uniquegenes(k,6),
            uniquegenes(k,6)= utrlocations(I(m),6);
        end;        
    end;    
end;


% Step 4:
% Find genes that come from one of the known chromosomes, and for which a valid UTR3 area exists.
%--------------------------------

I=find(isfinite(uniquegenes(:,2))&isfinite(uniquegenes(:,3)));
finalgenes=uniquegenes(I,:);
I=find(isnan(finalgenes(:,7))==0);
finalgenes=finalgenes(I,:);
