addpath ~/mlprojects/pol2rnaseq/matlab
cd ~/mlprojects/pol2rnaseq/matlab/
mex compute_pol2activityovergenes_c.c
cd ~/synergy_data/PolII/Mapping_results

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
n_chromosomes=length(chromosomenames);


%-----------------------------
% Read in empty region information
%-----------------------------

% read empty region chromosome ids and start and end locations
emptylocs=read_stringfile('~/synergy_data/analyses/empty_regions.txt',[32 9 10 13]);
n_emptyregions=length(emptylocs);


%-----------------------------
% Create bin array, one bin for each empty region
%-----------------------------
emptybininfo=zeros(length(emptylocs),6);
for i=1:length(emptylocs),
  emptytemp=emptylocs{i};

  % get chromosome index in chromosome list
  emptychromosome=emptytemp{1}(1:end);
  line_chrindex=-1;
  for k=1:n_chromosomes,
    if strcmp(emptychromosome,chromosomenames{k})==1,
      line_chrindex=k;
      break;
    end;
  end;

  emptybininfo(i,1)=line_chrindex; % chromosome index
  emptybininfo(i,2)=str2double(emptytemp{2}(1:end)); % start location
  emptybininfo(i,3)=str2double(emptytemp{3}(1:end)); % end location
  emptybininfo(i,4)=i;
  emptybininfo(i,5)=i;  % just a dummy ID for empty regions
  emptybininfo(i,6)=1;  % just a dummy strand ID for empty regions
end;



%-----------------------------
% Sort bin array in ascending order of start location. 
% The POL2 activity computation code assumes bins have been sorted like that.
%-----------------------------
[y,I1]=sort(emptybininfo(:,2));
emptybininfo=emptybininfo(I1,:);



%-----------------------------
% For each time point, count amount of POL2-read-basepairs that 
% overlap each bin, weighted by scores of the reads.
%-----------------------------

filenames={'pol0.mat','pol5.mat','pol10.mat','pol20.mat','pol40.mat','pol80.mat','pol160.mat','pol320.mat','pol640.mat','pol1280.mat'};
dvalues=[190 190 196 192 185 189 201 205 194 189];


% load up all data
allemptybins=cell(size(emptybininfo,1),10);
for timepoint=1:10,
  load(filenames{timepoint});
  d=dvalues(timepoint);
  max_duplicates=2;  
  subbin_length=200;
  
  for chr_index=1:n_chromosomes,
    timepoint
    chr_index
    I=find(emptybininfo(:,1)==chr_index);
    if length(I)>0,
      tempbins=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(emptybininfo(I,2)),int32(emptybininfo(I,3)),int8(emptybininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
      for k=1:length(I),
        allemptybins{I(k),timepoint}=tempbins{k};
      end;
    end;
  end;
  clear temppol;
end;

save all_empty_pol2bins.mat allemptybins emptybininfo


