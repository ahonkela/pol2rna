addpath ~/mlprojects/pol2rnaseq/matlab
cd ~/mlprojects/pol2rnaseq/matlab/
mex compute_pol2activityovergenes_c.c
cd ~/synergy_data/PolII/Mapping_results

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
n_chromosomes=length(chromosomenames);


%-----------------------------
% Read in gene information
%-----------------------------

% read gene ids and start and end locations
genelocs=read_stringfile('~/synergy_data/PolII/Mapping_results/genes.txt',[32 9 10 13]);
genelocs={genelocs{2:end}};
n_genes=length(genelocs);
geneids=zeros(n_genes,1);
for k=1:n_genes,
  geneids(k)=str2double(genelocs{k}{1}(6:end-1));
end;

% read gene ids and strand; match them to the previous gene list
genestrands_temp=read_stringfile('~/synergy_data/PolII/Mapping_results/genestrands_no_lrg_genes.txt',[32 9 10 13]);
geneids_temp=zeros(length(genestrands_temp),1);
for k=1:length(genestrands_temp),
  geneids_temp(k)=str2double(genestrands_temp{k}{3}(6:end-1));
end;
genestrands=zeros(n_genes,1);
for k=1:n_genes,
  l=find(geneids_temp==geneids(k));
  if length(l)==0,
    genestrands(k)=NaN;
  else    
    l=l(1);
    genestrands(k)=str2double(genestrands_temp{l}{2});
  end;
end;


%-----------------------------
% Create bin array, one bin for each gene
%-----------------------------
bininfo=zeros(length(genelocs),6);
for i=1:length(genelocs),
  genetemp=genelocs{i};

  % get chromosome index in chromosome list
  genechromosome=genetemp{2}(2:end-1);
  line_chrindex=-1;
  for k=1:n_chromosomes,
    if strcmp(genechromosome,chromosomenames{k})==1,
      line_chrindex=k;
      break;
    end;
  end;

  bininfo(i,1)=line_chrindex; % chromosome index
  bininfo(i,2)=str2double(genetemp{3}(2:end-1)); % start location
  bininfo(i,3)=str2double(genetemp{4}(2:end-1)); % end location
  bininfo(i,4)=i;
  bininfo(i,5)=geneids(i);
  bininfo(i,6)=genestrands(i);
end;



% Omit genes for which we can't find the strand or the chromosome...
I=find((isnan(bininfo(:,6))==0) | (bininfo(:,1)==-1));
bininfo=bininfo(I,:);



%-----------------------------
% Sort bin array in ascending order of start location. 
% The POL2 activity computation code assumes bins have been sorted like that.
%-----------------------------
[y,I1]=sort(bininfo(:,2));
bininfo=bininfo(I1,:);



%-----------------------------
% For each time point, count amount of POL2-read-basepairs that 
% overlap each bin, weighted by scores of the reads.
%-----------------------------


%allbins=zeros(n_genes,10);
filenames={'pol0.mat','pol5.mat','pol10.mat','pol20.mat','pol40.mat','pol80.mat','pol160.mat','pol320.mat','pol640.mat','pol1280.mat'};
dvalues=[190 190 196 192 185 189 201 205 194 189];


% load up all data
allgenebins=cell(size(bininfo,1),10);
for timepoint=1:10,
  load(filenames{timepoint});
  d=dvalues(timepoint);
  max_duplicates=2;  
  subbin_length=200;
  %binsthistime=zeros(n_genes,1);
  
  for chr_index=1:n_chromosomes,
    timepoint
    chr_index
    I=find(bininfo(:,1)==chr_index);
    tempbins=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int8(bininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
    for k=1:length(I),
      allgenebins{I(k),timepoint}=tempbins{k};
    end;
  end;
  clear temppol;
end;

save all_gene_pol2bins.mat allgenebins bininfo

binlengths=zeros(size(allgenebins,1),1);for k=1:size(allgenebins,1),binlengths(k)=length(allgenebins{k,1});end;





%-----------------------------------------------------------------
% Read in RNA activities (reads per kilobase of exon model)
%-----------------------------------------------------------------

filteredgenes=read_stringfile('rpkm_fdrfiltered.txt',[9 10 13 double(' ')]);
% Find POL indices corresponding to the filtered genes
filtered_genenames=cell(length(filteredgenes),1);
filtered_genenames_numeric=zeros(length(filteredgenes),1);

rna_summaryseries=nan*ones(size(bininfo,1),10);
for i=1:length(filteredgenes),
  gene_ensembleid=str2double(filteredgenes{i}{1}(6:end-1));
  bininfo_index=find(bininfo(:,5)==gene_ensembleid);
  if length(bininfo_index==1),
    for j=1:10,
      rna_summaryseries(bininfo_index,j)=str2double(filteredgenes{i}{j+1});
    end;  
  end;  
end;






%---------------------------------------------------
% normalize POL2 time points by overall read count.
% Do not normalize genes by their length, it does not
% make sense later when looking at sums over a fixed 
% number of bins.
%---------------------------------------------------
pol2overall_uniquehits = [101720334 91131727 72832980 93776268 93180719 94068955 87257588 85438329 92506713 94888353];
pol2filtered_hits = [82716552 76680010 51009272 78425211 74768277 73970821 59797272 56458506 83618025 81174490];
pol2filtered_multipliers = (pol2filtered_hits/mean(pol2filtered_hits)).^(-1);

allgenebins_normalized=cell(size(allgenebins,1),size(allgenebins,2));
for k=1:size(allgenebins,1),
  %genelength=bininfo(k,3)-bininfo(k,2)+1;
  genelength=1; % no gene length normalization...
  for l=1:size(allgenebins,2),
    allgenebins_normalized{k,l}=allgenebins{k,l}*pol2filtered_multipliers(l)/genelength;    
  end;
end;






%-----------------------------------------------------------------
% Create POL2 summary
%-----------------------------------------------------------------

% Version 1: sum of first few bins from the start site onwards

polsummarylength=20;
polsummaryoffset=5;

pol_summaryseries=nan*ones(size(bininfo,1),10);
nbinsinone=10;
for i=1:size(bininfo,1),

  % summarize by last bins 
  % (note that last bins = bins closest to transcription start)
  if (binlengths(i)>=polsummarylength+polsummaryoffset),
    for j=1:10,
      pol_summaryseries(i,j)=sum(allgenebins_normalized{i,j}((end-polsummarylength+1-polsummaryoffset):(end-polsummaryoffset)));
    end;
  end;
  
end;


% Version 2: sum of last few bins (towards the end of the gene)

polsummarylength=20;
polsummaryoffset=0;

pol_summaryseries=nan*ones(size(bininfo,1),10);
nbinsinone=10;
for i=1:size(bininfo,1),

  % summarize by first bins 
  % (note that first bins = bins closest to transcription end)
  if (binlengths(i)>=polsummarylength+polsummaryoffset),
    for j=1:10,
      pol_summaryseries(i,j)=sum(allgenebins_normalized{i,j}((1+polsummaryoffset):(polsummarylength+polsummaryoffset)));
    end;
  end;
  
end;







%---------------------------------------------------
% Find genes with: 
% - substantial POL2 activity 
% - substantial POL2 variance over time
% - substantial RNA activity 
% - substantial RNA variance over time
%---------------------------------------------------

pol_summaryseries_means=zeros(size(allgenebins_normalized,1),1);
for k=1:size(allgenebins_normalized,1),
  pol_summaryseries_means(k)=mean(pol_summaryseries(k,:));
end;
pol_summaryseries_variances=zeros(size(allgenebins_normalized,1),1);
for k=1:size(allgenebins_normalized,1),
  pol_summaryseries_variances(k)=var(pol_summaryseries(k,:));
end;
rna_summaryseries_means=zeros(size(allgenebins_normalized,1),1);
for k=1:size(allgenebins_normalized,1),
  rna_summaryseries_means(k)=mean(rna_summaryseries(k,:));
end;
rna_summaryseries_variances=zeros(size(allgenebins_normalized,1),1);
for k=1:size(allgenebins_normalized,1),
  rna_summaryseries_variances(k)=var(rna_summaryseries(k,:));
end;

Isubstantialpol2=find((pol_summaryseries_means>=exp(9)) & ...
		      (pol_summaryseries_variances>=exp(15)) & ...
		      (rna_summaryseries_means>=exp(1.3)) & ...
		      (rna_summaryseries_variances>=exp(1.7)));



%---------------------------------------------------
% Create a clustering of the joint POL2-RNA series, restricted to
% genes with substantial activity
%---------------------------------------------------

%jointseries=[pol_summaryseries(Isubstantialpol2,:) rna_summaryseries(Isubstantialpol2,:)];

jointseries=[pol_summaryseries(Isubstantialpol2,:) rna_summaryseries_derivatives(Isubstantialpol2,:)];

% normalize over time...
for i=1:size(jointseries,1),
  jointseries(i,1:10)=jointseries(i,1:10)/sum(jointseries(i,1:10));
  jointseries(i,11:20)=jointseries(i,11:20)/sum(abs(jointseries(i,11:20)));
end;
nclusters=20;
idx=kmeans(jointseries,nclusters,'EmptyAction','singleton','MaxIter',1000);

nmembers=zeros(nclusters,1);
for i=1:nclusters,
  I=find(idx==i);
  nmembers(i)=length(I);
end;
[y,I2]=sort(-nmembers);

for k=1:5,
  figure;
  for j0=1:4,
    subplot(2,2,j0); j=(k-1)*4+j0;

    I=find(idx==I2(j));

    clustmean_pol=mean(jointseries(I,1:10));
    clustvar_pol=var(jointseries(I,1:10));
    clustmean_rna=mean(jointseries(I,11:20));
    clustvar_rna=var(jointseries(I,11:20));
    
    boxplot(jointseries(I,1:10),'colors',[1 0.5 0.5]);
    hold on;
    boxplot(jointseries(I,11:20),'colors',[0.5 0.5 1]);
    hold on; 
    h=plot(clustmean_pol,'r-');
    set(h,'LineWidth',2);
    h=plot(clustmean_rna,'b-');
    set(h,'LineWidth',2);
    
    mytitle=sprintf('POL2(red)-RNA(blue) cluster %d (%d members)', j, length(I));
    title(mytitle);
  end;
end;






delaytries=[0:0.5:40];
meancorrelations=zeros(length(delaytries),1);
ncandidates=zeros(length(delaytries),1);

for delayindex=1:length(delaytries),

%-----------------------------------------------------------------
% Create simple estimate of RNA derivative, by fitting splines.
% Apply a delay, for computation of the correlations.
%-----------------------------------------------------------------

timepoints=[0 5 10 20 40 80 160 320 640 1280];
rna_delay=delaytries(delayindex)
spline_spacing=0.1;

rna_summaryseries_derivatives=nan*ones(size(bininfo,1),10);
for i=1:size(bininfo,1),
  if sum(isnan(rna_summaryseries(i,:)))==0,
    tempseries=rna_summaryseries(i,:);
    tempseries=tempseries/sum(tempseries);
    interpolatedpoints=[0:spline_spacing:(1281+ceil(rna_delay))];
    firstderivative=0;
    lastderivative=0;
    % lastderivative=(tempseries(end)-tempseries(end-1))/(1280-640);
    tempspline=spline(timepoints,[firstderivative tempseries lastderivative],interpolatedpoints);
    for j=1:10,
      tempI=find(round((timepoints(j)+rna_delay)/spline_spacing)*spline_spacing==interpolatedpoints);
      rna_summaryseries_derivatives(i,j)=(tempspline(tempI(1)+1)-tempspline(tempI(1)))/spline_spacing;
    end;  
  end;
end;






%---------------------------------------------------
% Compute correlations between POL2 and (delayed) RNA derivatives
%---------------------------------------------------
lengthconsidered=7;
simpledelay=0;

corrs=zeros(length(Isubstantialpol2),1);
pvals=zeros(length(Isubstantialpol2),1);
for i=1:length(corrs),
  if (var(pol_summaryseries(Isubstantialpol2(i),1:lengthconsidered))>0) && (var(rna_summaryseries_derivatives(Isubstantialpol2(i),(1+simpledelay):(lengthconsidered+simpledelay)))>0),
    [tempcorr,pval]=corr([pol_summaryseries(Isubstantialpol2(i),1:lengthconsidered)' rna_summaryseries_derivatives(Isubstantialpol2(i),(1+simpledelay):(lengthconsidered+simpledelay))']);
  else
    tempcorr=[0 0;0 0];
    pval=[1 1;1 1];
  end;  
  corrs(i)=tempcorr(1,2);
  pvals(i)=pval(1,2);
end;

[y,I3]=sort(-corrs);
candidates=find((corrs>0).*(pvals <= 0.01));

meancorrelations(delayindex)=mean(corrs);
ncandidates(delayindex)=length(candidates);


end;



figure; 
subplot(2,1,1); plot(delaytries,meancorrelations); 
xlabel('RNA delay (min)'); ylabel('Mean POL-RNA corr'); 
title('RNAderivative-POL correlations, POL=sum of [0bp-4000bp] from gene end')
subplot(2,1,2); plot(delaytries,ncandidates); 
xlabel('RNA delay (min)'); ylabel('N.of significant pos.corrs.');





%timepoint=4;
%chr_index=2;
%  load(filenames{timepoint});
%  d=dvalues(timepoint);
%  max_duplicates=2;  
%  subbin_length=200;
%  binsthistime=zeros(n_genes,1);
%    I=find(bininfo(:,1)==chr_index);
%    
%cd ~/mlprojects/pol2rnaseq/matlab/
%mex compute_pol2activityovergenes_c.c
%cd ~/synergy_data/PolII/Mapping_results
%
%    tempbins=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int8(bininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));


    
    
%-----------------------------------------------------------------
% compute a "spatial correlation matrix with delays" in a really 
% simplistic way, using simple interpolation for delays that do
% not fall exactly onto a particular time point
% assumption: x(end) at time t correlates with 
%             x(end-k) at time t-delay*k
%-----------------------------------------------------------------    

binlengths=zeros(size(allgenebins,1),1);for k=1:size(allgenebins,1),binlengths(k)=length(allgenebins{k,1});end;

profilelength=100;
%I2=find(binlengths>=profilelength);
I2=find((binlengths>=profilelength) & (pol_sumseries_means>=10) & (pol_sumseries_variances>=10));
%I2=[6586]; % GREB1

%tempallgenebins=allgenebins_normalized(I2,:);
tempallgenebins=allgenebins_normalized(I2,:);
for k=1:size(tempallgenebins,1),
  for l=1:size(tempallgenebins,2),
    tempallgenebins{k,l}=tempallgenebins{k,l}(end-profilelength+1:end);
  end;
end;
% aggregate bins, 10 bins in one
nbinsinone=10;
aggregatedlength=profilelength/nbinsinone;
for k=1:size(tempallgenebins,1),
  for l=1:size(tempallgenebins,2),
    tempbins=zeros(aggregatedlength,1);
    for m=1:aggregatedlength,
      tempbins(m)=sum(tempallgenebins{k,l}((1+(m-1)*nbinsinone):(nbinsinone+(m-1)*nbinsinone)));
    end;    
    tempallgenebins{k,l}=tempbins;
  end;
end;


addpath ~/mlprojects/pol2rnaseq/matlab/
cd ~/mlprojects/pol2rnaseq/matlab/
mex compute_fakecorrelationmatrix_c.c
cd ~/synergy_data/PolII/Mapping_results

speeds=[[-4:0.1:4]];
entropies=zeros(length(speeds),1);
sumpearsons=zeros(length(speeds),1);
for k=1:length(speeds),
  pol_spatialspeed=speeds(k); %(2000/5)/200;

  % compute correlation matrix
  [binmeans,binvariances,binsamples,corrmatrix,corrsamples]=compute_fakecorrelationmatrix_c(int32(size(tempallgenebins,1)),int32(10),int32(aggregatedlength),tempallgenebins,double(pol_spatialspeed));

  tempmatrix=corrmatrix-(binmeans*binmeans');
  entropies(k)=0.5*aggregatedlength*log(2*pi*exp(1))+0.5*(log(det(tempmatrix(1:aggregatedlength,1:aggregatedlength)/tempmatrix(1,1)))+log(tempmatrix(1,1)));

  tempmatrix2=(corrmatrix-(binmeans*binmeans'))./sqrt(binvariances*binvariances');
  sumpearsons(k)=sum(sum(abs(tempmatrix2(1:aggregatedlength,1:aggregatedlength))));
  
%  figure; 
%  %imagesc((corrmatrix-(binmeans*binmeans'))./sqrt(binvariances*binvariances')); 
%  imagesc((corrmatrix-(binmeans*binmeans'))); 
%  colormap(gray);
%  title(sprintf('lag (min/200bp) %f, entropy %f',pol_spatialspeed, entropies(k)));  
end;




















%-----------------------------------------------------------------
% compute a mean normalized gene profile, weighted by sum of bins
%-----------------------------------------------------------------
profilelength=100;
I2=find(binlengths>=profilelength);



profile=zeros(length(I2),profilelength);
for k=1:length(I2),
  profile(k,:)=tempbins{I2(k)}(end:-1:end-profilelength+1);
end;

normalizedprofile=[];
profilemean=zeros(1,profilelength);
normalizer=0;
nsamples=0;
for k=1:length(I2),
  if sum(profile(k,:))>0,
    tempprofile=profile(k,:)/sum(profile(k,:));
    normalizedprofile=[normalizedprofile;tempprofile];
    profilemean=profilemean+tempprofile*sum(profile(k,:));
    normalizer=normalizer+sum(profile(k,:));
    nsamples=nsamples+1;
  end;  
end;
profilemean=profilemean/normalizer;


profilevar=zeros(1,profilelength);
for k=1:length(I2),
  if sum(profile(k,:))>0,
    profilevar=profilevar+ sum(profile(k,:))*( ( profile(k,:)/sum(profile(k,:)) ) -profilemean).^2;  
  end;    
end;
profilevar=profilevar/normalizer;    

clf;
boxplot(normalizedprofile); hold on;
h=plot(profilemean,'g-'); set(h,'LineWidth',2);
h=plot(profilemean-sqrt(profilevar),'c-'); set(h,'LineWidth',2);
h=plot(profilemean+sqrt(profilevar),'c-'); set(h,'LineWidth',2);
title(sprintf('Normalized Pol2-profile from transcription star site\n(Weighted mean and STD over 1177 genes, weighted by their Pol2 activity)'))
xlabel('Bins from gene start, bins size=200 basepairs');
ylabel('Normalized Pol2-activity at each bin');









################### OLD CODE ########################

% compute bin means
binmeans=zeros(profilelength+1,1);
binsamples=zeros(profilelength+1,1);
for i=1:profilelength, % omit RNA for now
  i
  for g=1:length(I2),
    for timeindex=1:10,
      t=mytimepoints(timeindex);
      value_t=allgenebins{I2(g),timeindex}(i);
      binmeans(i)=binmeans(i)+value_t;
      binsamples(i)=binsamples(i)+1;
    end;
  end;  
end;
binmeans=binmeans./binsamples;

% compute bin variances
binvariances=zeros(profilelength+1,1);
for i=1:profilelength, % omit RNA for now
  i
  for g=1:length(I2),
    for timeindex=1:10,
      t=mytimepoints(timeindex);
      value_t=allgenebins{I2(g),timeindex}(i);
      binvariances(i)=binvariances(i)+(value_t-binmeans(i))^2;
    end;
  end;  
end;
binvariances=binvariances./binsamples;

% spatial speed of POL2 (with pauses and all) in bins/minute = (basepairs/minute)/binlength
pol_spatialspeed=(2000/5)/200;

corrmatrix=zeros(profilelength+1,profilelength+1);
ncorrsamples=zeros(profilelength+1,profilelength+1);
mytimepoints=[0 5 10 20 40 80 160 320 640 1280];
for timeindex=1:10,
  t=mytimepoints(timeindex);
  for g=1:length(I2),
    [t g]
    for i=1:profilelength, % omit RNA for now
      value_t=allgenebins{I2(g),timeindex}(i);
            
      for j=1:profilelength, % omit RNA for now
	t2=t+(j-i)/pol_spatialspeed;
	if (t2 >= 0) & (t2 <= 1280),

	  timeindex2=max(find(mytimepoints<=t2));
	  if (timeindex2<10),
	    value_t2=allgenebins{I2(g),timeindex2}(j)+(t2-mytimepoints(timeindex2))/(mytimepoints(timeindex2+1)-mytimepoints(timeindex2))*(allgenebins{I2(g),timeindex2+1}(j)-allgenebins{I2(g),timeindex2}(j));
	  else
	    value_t2=allgenebins{I2(g),timeindex2}(j);
	  end;
	  
	  corrmatrix(i,j)=corrmatrix(i,j)+value_t*value_t2;
	  ncorrsamples(i,j)=ncorrsamples(i,j)+1;
	end;
      end;      
    end;
    
  end;
end;











###################### debug code


allbins=zeros(n_genes,10);
filenames={'pol0.mat','pol5.mat','pol10.mat','pol20.mat','pol40.mat','pol80.mat','pol160.mat','pol320.mat','pol640.mat','pol1280.mat'};
dvalues=[190 190 196 192 185 189 201 205 194 189];


for timepoint=1:10,
  load(filenames{timepoint});
  d=dvalues(timepoint);
  max_duplicates=2;  
  binsthistime=zeros(n_genes,1);
  
  for chr_index=2:2,
    timepoint
    chr_index
    I=6602;
    tempbins=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
    binsthistime(I)=tempbins;
  end;
  
  allbins(:,timepoint)=binsthistime;
  clear temppol;
end;






######################## old code ########################


load pol5.mat
d=190;
max_duplicates=2;
bins5=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
clear temppol;

load pol10.mat
d=196;
max_duplicates=2;
bins10=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
clear temppol;

load pol20.mat
d=192;
max_duplicates=2;
bins20=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
clear temppol;

load pol40.mat
d=185;
max_duplicates=2;
bins40=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
clear temppol;

load pol80.mat
d=189;
max_duplicates=2;
bins80=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
clear temppol;

load pol160.mat
d=201;
max_duplicates=2;
bins160=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
clear temppol;

load pol320.mat
d=205;
max_duplicates=2;
bins320=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
clear temppol;

load pol640.mat
d=194;
max_duplicates=2;
bins640=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
clear temppol;

load pol1280.mat
d=189;
max_duplicates=2;
bins1280=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
clear temppol;

bins=[bins0 bins5 bins10 bins20 bins40 bins80 bins160 bins320 bins640 bins1280];
save /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/pol2bins.mat bins


%-----------------------------
% Code for converting .bed files of POL2 data into matlab cell
% arrays. Only needs to be run once.
%-----------------------------

if 0,
!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_No_Treat_PolII_unique.bed.gz
temppol=read_mappingfile('MCF7_No_Treat_PolII_unique.bed');
save pol0.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_No_Treat_PolII_unique.bed
temppol=1;

!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_5min_PolII_unique.bed.gz
temppol=read_mappingfile('MCF7_E2_5min_PolII_unique.bed');
save pol5.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_5min_PolII_unique.bed
temppol=1;

!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_10min_PolII_rep_unique.bed.gz
temppol=read_mappingfile('MCF7_E2_10min_PolII_rep_unique.bed');
save pol10.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_10min_PolII_rep_unique.bed
temppol=1;

!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_20min_PolII_unique.bed.gz
temppol=read_mappingfile('MCF7_E2_20min_PolII_unique.bed');
save pol20.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_20min_PolII_unique.bed
temppol=1;

!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_40min_PolII_unique.bed.gz
temppol=read_mappingfile('MCF7_E2_40min_PolII_unique.bed');
save pol40.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_40min_PolII_unique.bed
temppol=1;

!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_80min_PolII_unique.bed.gz
temppol=read_mappingfile('MCF7_E2_80min_PolII_unique.bed');
save pol80.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_80min_PolII_unique.bed
temppol=1;

!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_160min_PolII_unique.bed.gz
temppol=read_mappingfile('MCF7_E2_160min_PolII_unique.bed');
save pol160.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_160min_PolII_unique.bed
temppol=1;

!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_320min_PolII_unique.bed.gz
temppol=read_mappingfile('MCF7_E2_320min_PolII_unique.bed');
save pol320.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_320min_PolII_unique.bed
temppol=1;

!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_640min_PolII_unique.bed.gz
temppol=read_mappingfile('MCF7_E2_640min_PolII_unique.bed');
save pol640.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_640min_PolII_unique.bed
temppol=1;

!/bin/gunzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_1280min_PolII_unique.bed.gz
temppol=read_mappingfile('MCF7_E2_1280min_PolII_unique.bed');
save pol1280.mat temppol
!/bin/gzip /home/jaakkopeltonen/synergy_data/PolII/Mapping_results/MCF7_E2_1280min_PolII_unique.bed
temppol=1;

end;
