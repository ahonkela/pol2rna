addpath ~/mlprojects/pol2rnaseq/matlab
cd ~/mlprojects/pol2rnaseq/matlab/
mex compute_pol2activityovergenes_c.c
cd ~/synergy_data/PolII/Mapping_results

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
n_chromosomes=length(chromosomenames);

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
% The activity computation code assumes bins have been sorted like that.
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

    
    
%-----------------------------------------------------------------
% compute a "spatial correlation matrix with delays" in a really 
% simplistic way, using simple interpolation for delays that do
% not fall exactly onto a particular time point
% assumption: x(end) at time t correlates with 
%             x(end-k) at time t-delay*k
%-----------------------------------------------------------------    


% compute correlation matrix
profilelength=100;
I2=find(binlengths>=profilelength);



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
