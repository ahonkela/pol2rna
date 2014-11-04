addpath /share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/pol2rnaseq/matlab

load /share/mi/workspace/jtpelto/synergy/synergy_data/analyses/new_bininfo/bininfo_dec2012.mat

A=load('transcript_exonpos.mat');
A=A.A;
tcountsfile=read_stringfile('active_transcripts.txt',[' ' ',' 9 10 13]); % space comma tab linefeed carriagereturn
tcounts=zeros(length(tcountsfile),length(tcountsfile{1}));
for i=1:length(tcountsfile),
  tcounts(i,1)=str2double(tcountsfile{i}{1}(5:end)); % gene id
  tcounts(i,2)=str2double(tcountsfile{i}{2}(5:end)); % transcript id
  tcounts(i,3)=str2double(tcountsfile{i}{3}); % mean number of reads
end;

% Find all genes where there are active transcripts
geneids=unique(tcounts(:,1));

% For all those genes, find the active transcripts and
% their corresponding minimum and maximum basepair position
genestarts=zeros(length(geneids),1);
geneends=zeros(length(geneids),1);
nproblemtranscripts=0;
problemtranscripts=[];
% Istartexons=find(A(:,5)==1);
% Astartexons=A(Istartexons,:);
for k=1:length(geneids),
  if mod(k,1000)==0,
    k
  end;
  gid=geneids(k);
  Itrids=find(tcounts(:,1)==gid);
  trids=tcounts(Itrids,2);
  % Find start exon positions of each transcript
  Iexons=find(ismember(A(:,2),trids));
  if length(Iexons)>0,
    genestarts(k)=min(A(Iexons,3));
    geneends(k)=max(A(Iexons,4));
  else
    genestarts(k)=nan;
    geneends(k)=nan;
  end;
  if length(Iexons)<length(trids),
    missingtrs=setdiff(trids,A(Iexons,2));
    nproblemtranscripts=nproblemtranscripts+length(missingtrs);
    problemtranscripts=[problemtranscripts; missingtrs];
  end;
end;

genes_without_active_uptodate_transcripts=geneids(find(isnan(genestarts)==1));

