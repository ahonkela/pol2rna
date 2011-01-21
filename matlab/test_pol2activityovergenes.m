chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrMT'};
n_chromosomes=length(chromosomenames);

genelocs=read_stringfile('~/synergy_data/PolII/Mapping_results/genes.txt',[32 9 10 13]);
genelocs={genelocs{2:end}};
n_genes=length(genelocs);


%-----------------------------
% Create bin array, one bin for each gene
%-----------------------------

bininfo=zeros(length(genelocs),4);
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
end;

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


allbins=zeros(n_genes,10);
filenames={'pol0.mat','pol5.mat','pol10.mat','pol20.mat','pol40.mat','pol80.mat','pol160.mat','pol320.mat','pol640.mat','pol1280.mat'};
dvalues=[190 190 196 192 185 189 201 205 194 189];


for timepoint=1:10,
  load(filenames{timepoint});
  d=dvalues(timepoint);
  max_duplicates=2;  
  binsthistime=zeros(n_genes,1);
  
  for chr_index=1:n_chromosomes,
    timepoint
    chr_index
    I=find(bininfo(:,1)==chr_index);
    tempbins=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int32(d),int32(max_duplicates));
    binsthistime(I)=tempbins;
  end;
  
  allbins(:,timepoint)=binsthistime;
  clear temppol;
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
