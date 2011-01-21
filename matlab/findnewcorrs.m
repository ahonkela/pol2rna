%geneseries=zeros(length(genelocs),7);
%for i=1:length(genelocs),
% for k=1:7,
% geneseries(i,k)=str2double(genelocs{i}{4+k}(2:end-1));
%end;
%end;
geneseries=double(rawcounts);

%chr_index=1;
%Ichr=find(bininfo(:,1)==chr_index);
%geneindices=bininfo(Ichr,4);
%genebins=geneseries(geneindices,:);


pol2overall_uniquehits = [101720334 91131727 72832980 93776268 93180719 94068955 87257588 85438329 92506713 94888353];
pol2filtered_hits = [82716552 76680010 51009272 78425211 74768277 73970821 59797272 56458506 83618025 81174490];
pol2filtered_multipliers = pol2filtered_hits/mean(pol2filtered_hits);
polbins_normalized=allbins;
for i=1:10,
  polbins_normalized(:,i)=allbins(:,i)*pol2filtered_multipliers(i);
end;


rna_uniquehits = [21695936 20435564 22212154 21769755 23547458 ...
		  22137760 20776819 26753863 25437235 26747370];
rna_multipliers = rna_uniquehits/mean(rna_uniquehits);
rnabins_normalized=geneseries;
for i=1:7,
  rnabins_normalized(:,i)=geneseries(:,i)*rna_multipliers(i);
end;

% sort in the same ascending order of start location as in the pol-bins
rnabins_normalized=rnabins_normalized(I1,:);



 

%I5=find([(sum((rnabins_normalized),2)>51).* ...
%	 (sum((rnabins_normalized==0),2)==0).*(sum((polbins_normalized==0),2)==0)]);

I5=find([(sum((rnabins_normalized==0),2)==0).*(sum((polbins_normalized==0),2)==0)]);
rnabins_nonzero=rnabins_normalized(I5,:);
polbins_nonzero=polbins_normalized(I5,:);

geneindices_nonzeropairs=I1(I5);



lengthconsidered=4;
delay=1;
corrs=zeros(length(I5),1);
pvals=zeros(length(I5),1);
for i=1:length(corrs),
  [tempcorr,pval]=corr([polbins_nonzero(i,1:lengthconsidered)' rnabins_nonzero(i,(1+delay):(lengthconsidered+delay))']);
  corrs(i)=tempcorr(1,2);
  pvals(i)=pval(1,2);
end;
[y,I3]=sort(-corrs);
candidates=find((corrs>0).*(pvals <= 0.01));
candidate_geneindices = geneindices_nonzeropairs(candidates);



corrs=zeros(length(I5),1);
pvals=zeros(length(I5),1);
for i=1:length(corrs),
  [tempcorr,pval]=corr([polbins_nonzero(i,1:4)' rnabins_nonzero(i,1:4)']);
  corrs(i)=tempcorr(1,2);
  pvals(i)=pval(1,2);
end;
[y,I3]=sort(-corrs);
candidates=find((corrs>0).*(pvals <= 0.01));
candidate_geneindices = geneindices_nonzeropairs(candidates);




% significance testing
ntrials=10000;
trialresults=zeros(ntrials,1);

%[allcorrs,allpvals]=corr(polbins_nonzero', rnabins_nonzero');  

for trial=1:ntrials,
  if mod(trial,10)==0, 
    trial
  end;
  myperm=randperm(length(I5));
  tempcorrs=zeros(length(I5),1);
  temppvals=zeros(length(I5),1);   
  for i=1:length(tempcorrs),
    [tempcorr,temppval]=corr(polbins_nonzero(i,1:10)', rnabins_nonzero(myperm(i),1:10)');    
    tempcorrs(i)=tempcorr;
    temppvals(i)=temppval;
  end;
  tempcandidates=find((tempcorrs>0).*(temppvals <= 0.01));
  trialresults(trial)=length(tempcandidates);
end;




for k=1:length(genelocs),
  geneids(k)=str2double(genelocs{k}{1}(7:end-1));
end;
% 17578 = GREB1

