load all_gene_pol2bins.mat



%filteredgenes=read_stringfile('filteredgenes.txt',[9 10 13 double(' ')]);
filteredgenes=read_stringfile('fdrfiltereddata.txt',[9 10 13 double(' ')]);
% Find POL indices corresponding to the filtered genes
filtered_genenames=cell(length(filteredgenes),1);
filtered_genenames_numeric=zeros(length(filteredgenes),1);
for i=1:length(filteredgenes),
  filtered_genenames{i}=filteredgenes{i}{1}(2:end-1);
  filtered_genenames_numeric(i)=str2double(filtered_genenames{i}(5:end));
end;
filtered_geneseries=zeros(length(filteredgenes),10);
for i=1:length(filteredgenes),
  for j=1:10,
    filtered_geneseries(i,j)=str2double(filteredgenes{i}{j+1});
  end;  
end;

allrnabins=nan*ones(size(bininfo,1),10);
for i=1:length(filteredgenes),    
  i
  % find corresponding gene index in the pol array
  polindex=find(bininfo(:,5)==filtered_genenames_numeric(i));
  if size(polindex,1)==1,
    allrnabins(polindex,:)=filtered_geneseries(i,:);
end;



%temp_polnames_numeric=zeros(size(allbins,1),1);
%for i=1:length(I1),
%  temp_polnames_numeric(i)=str2double(genelocs{I1(i)}{1}(6:(end-1)));
%end;

filtered_polseries=-1*ones(length(filteredgenes),10);




%pol_gene
%
%geneseries=double(rawcounts);

%geneseries=zeros(length(genelocs),7);
%for i=1:length(genelocs),
% for k=1:7,
% geneseries(i,k)=str2double(genelocs{i}{4+k}(2:end-1));
%end;
%end;

%chr_index=1;
%Ichr=find(bininfo(:,1)==chr_index);
%geneindices=bininfo(Ichr,4);
%genebins=geneseries(geneindices,:);



pol2overall_uniquehits = [101720334 91131727 72832980 93776268 93180719 94068955 87257588 85438329 92506713 94888353];
pol2filtered_hits = [82716552 76680010 51009272 78425211 74768277 73970821 59797272 56458506 83618025 81174490];
pol2filtered_multipliers = pol2filtered_hits/mean(pol2filtered_hits);
polbins_normalized=zeros(length(filtered_genenames),10);
for i=1:10,
  polbins_normalized(:,i)=filtered_polseries(:,i)*pol2filtered_multipliers(i);
end;


rna_uniquehits = [21695936 20435564 22212154 21769755 23547458 ...
		  22137760 20776819 26753863 25437235 26747370];
rna_multipliers = rna_uniquehits/mean(rna_uniquehits);
rnabins_normalized=zeros(length(filtered_genenames),10);
for i=1:10,
  rnabins_normalized(:,i)=filtered_geneseries(:,i)*rna_multipliers(i);
end;

% compute variance
rnanormalized_variance=zeros(length(filtered_genenames),1);
for i=1:length(filtered_genenames),
  rnanormalized_variance(i)=var(rnabins_normalized(i,:));
end;

% compute variance
polnormalized_variance=zeros(length(filtered_genenames),1);
for i=1:length(filtered_genenames),
  polnormalized_variance(i)=var(polbins_normalized(i,:));
end;


% sort in the same ascending order of start location as in the pol-bins
% rnabins_normalized=rnabins_normalized(I1,:);




[N,X]=hist(log(sum(rnabins_normalized,2)+1),60);
semilogx(exp(X)-1,N);
xlabel('sum of RNA (read hits) over all time points')
ylabel('number of genes with this sum');
title('RNA, histogram of overall sums');

[N,X]=hist(log(sum(polbins_normalized,2)+1),60);
semilogx(exp(X)-1,N);
xlabel('sum of POL2 (basepair hits) over all time points')
ylabel('number of genes with this sum');
title('POL2, histogram of overall sums');

[N,X]=hist(log(rnanormalized_variance),60);
semilogx(exp(X),N);
xlabel('variance of RNA time series')
ylabel('number of genes with this variance');
title('RNA, variance of time series');

[N,X]=hist(log(polnormalized_variance+1),60);
semilogx(exp(X)-1,N);
xlabel('variance of POL time series')
ylabel('number of genes with this variance');
title('POL2, variance of time series');


loglog(sum(polbins_normalized,2),sum(rnabins_normalized,2),'k.')
xlabel('sum of POL2 (basepair hits) over all time points')
ylabel('sum of RNA (read hits) over all time points')

 

%I5=find([(sum((rnabins_normalized),2)>51).* ...
%	 (sum((rnabins_normalized==0),2)==0).*(sum((polbins_normalized==0),2)==0)]);

%I5=find([(sum((rnabins_normalized==0),2)==0).*(sum((polbins_normalized==0),2)==0)]);
%I5=find([(sum((rnabins_normalized==0),2)==0).*(sum((polbins_normalized==0),2)==0)]);




I5=find((sum(polbins_normalized,2)>1e5).*(rnanormalized_variance>1e2).*(polnormalized_variance>1e5));
rnabins_nonzero=rnabins_normalized(I5,:);
polbins_nonzero=polbins_normalized(I5,:);

%geneindices_nonzeropairs=I1(I5);

lengthconsidered=7;
delay=1;

corrs=zeros(length(I5),1);
pvals=zeros(length(I5),1);
for i=1:length(corrs),
  if (var(polbins_nonzero(i,1:lengthconsidered))>0) && (var(rnabins_nonzero(i,(1+delay):(lengthconsidered+delay)))>0),
    [tempcorr,pval]=corr([polbins_nonzero(i,1:lengthconsidered)' rnabins_nonzero(i,(1+delay):(lengthconsidered+delay))']);
  else
    tempcorr=[0 0;0 0];
    pval=[1 1;1 1];
  end;  
  corrs(i)=tempcorr(1,2);
  pvals(i)=pval(1,2);
end;

[y,I3]=sort(-corrs);
candidates=find((corrs>0).*(pvals <= 0.01));
%candidate_geneindices = geneindices_nonzeropairs(candidates);



%corrs=zeros(length(I5),1);
%pvals=zeros(length(I5),1);
%for i=1:length(corrs),
%  [tempcorr,pval]=corr([polbins_nonzero(i,1:4)' rnabins_nonzero(i,1:4)']);
%  corrs(i)=tempcorr(1,2);
%  pvals(i)=pval(1,2);
%end;
%[y,I3]=sort(-corrs);
%candidates=find((corrs>0).*(pvals <= 0.01));
%candidate_geneindices = geneindices_nonzeropairs(candidates);




% significance testing
%[allcorrs,allpvals]=corr(polbins_nonzero', rnabins_nonzero');  

ntrials=100;
trialresults=zeros(ntrials,1);

for trial=1:ntrials,
  if mod(trial,10)==0, 
    trial
  end;
  myperm=randperm(length(I5));
  tempcorrs=zeros(length(I5),1);
  temppvals=zeros(length(I5),1);   
  for i=1:length(tempcorrs),
    if (var(polbins_nonzero(i,1:lengthconsidered))>0) && (var(rnabins_nonzero(myperm(i),(1+delay):(lengthconsidered+delay)))>0),
      [tempcorr,temppval]=corr(polbins_nonzero(i,1:lengthconsidered)', rnabins_nonzero(myperm(i),(1+delay):(lengthconsidered+delay))');    
    else
      tempcorr=0;
      temppval=1;
    end;
    
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

