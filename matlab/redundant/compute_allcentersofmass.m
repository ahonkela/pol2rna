analysisdir='/Users/hasanogul/jaakkos_files/synergy/analyses/';
stringfile=read_stringfile([analysisdir 'earlygenes_t5.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t5=tempgenes;
stringfile=read_stringfile([analysisdir 'earlygenes_t10.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t10=tempgenes;
stringfile=read_stringfile([analysisdir 'earlygenes_t20.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t20=tempgenes;
stringfile=read_stringfile([analysisdir 'earlygenes_t40.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t40=tempgenes;
stringfile=read_stringfile([analysisdir 'earlygenes_t80.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t80=tempgenes;
interestinggenes=unique([interestinggenes_t5;interestinggenes_t10;interestinggenes_t20;interestinggenes_t40;interestinggenes_t80]);


masscenters1=zeros(length(interestinggenes),10);
masscenters2=zeros(length(interestinggenes),10);

for m=1:length(interestinggenes),
m
  ensemblid=interestinggenes(m);

  gene_index=find(bininfo(:,5)==ensemblid);
  if isempty(gene_index),
    masscenters1(m,:)=nan;
    masscenters2(m,:)=nan;
    continue;
  end;

  if length(allgenebins{gene_index,1})<20,
    masscenters1(m,:)=nan;
    masscenters2(m,:)=nan;
    continue;
  end;

  tempprofile=[];
  for k=1:10,
    tempprofile(k,:)=allgenebins{gene_index,k}(end:-1:1);
  end;
  if min(sum(tempprofile,2))==0,
    masscenters1(m,:)=nan;
    masscenters2(m,:)=nan;
    continue;
  end;


  % apply normalization of time points (geommeanmedian, based on last time points)
  normalized_tempprofile1=[];
  for k=1:10,
    normalized_tempprofile1(k,:)=tempprofile(k,:)*normalizernumbers(k);
  end;

  % compute mean profile
  meanprofile=mean(normalized_tempprofile1,1);

  % apply normalization by mean profile
  normalized_tempprofile2=[];
  I = find(meanprofile>0);
  for k=1:10,
    normalized_tempprofile2(k,:)=normalized_tempprofile1(k,:);
    normalized_tempprofile2(k,I)=normalized_tempprofile1(k,I)./meanprofile(I);
  end;

  % compute center of mass, first without normalization with the mean profile
  centerofmass1=[];
  for k=1:10,
    centerofmass1(k)=sum(normalized_tempprofile1(k,:).*[1:size(normalized_tempprofile1,2)]*200)/...
      sum(normalized_tempprofile1(k,:));
  end;
  masscenters1(m,:)=centerofmass1;

  % compute center of mass, using normalization with the mean profile
  centerofmass2=[];
  for k=1:10,
    centerofmass2(k)=sum(normalized_tempprofile2(k,:).*[1:size(normalized_tempprofile2,2)]*200)/...
      sum(normalized_tempprofile2(k,:));
  end;
  masscenters2(m,:)=centerofmass2;
end;
