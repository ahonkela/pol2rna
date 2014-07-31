mybasedir='/Users/hasanogul/jaakkos_files/synergy/mlprojects/';
pol2dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/PolII/Mapping_results/'

% for kernel-level computations
path1=[mybasedir 'kern/matlab/jaakko_testversion/']
% for model-level computations
path2=[mybasedir 'gpsim/matlab/jaakko_testversion/'];
% for optimiDefaultConstraint.m
path3=[mybasedir 'optimi/matlab/'];
% for lnDiffErfs.m
path4=[mybasedir 'ndlutil/matlab/'];
% for addPrior.m
path5=[mybasedir 'prior/matlab/'];
% for dist2.m
path6=[mybasedir 'matlab/netlab/NETLAB3p3/'];
% for modelTieParam.m
path7=[mybasedir 'mltools/matlab/'];
% for plotting etc. functions
path8=[mybasedir 'pol2rnaseq/matlab/'];

addpath(path1,path2,path3,path4,path5,path6,path7,path8)

cd(pol2dir)
load pol2_for_matti_ver2.mat


pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_last_0kb_to_4kb;
pol_summaryseries_normalized=normalizedbyfilteredcount_pol_summaryseries_last_0kb_to_4kb;
normalizernumbers=pol_summaryseries_normalized./pol_summaryseries_unnormalized;

% just take the normalizers from one example gene, they're the
% same for all genes up to computational precision
normalizernumbers=normalizernumbers(2000,:); 

load all_gene_pol2bins.mat

ensemblid=196208; % GREB1
gene_index=find(bininfo(:,5)==ensemblid);
tempprofile=[];
for k=1:10,
  tempprofile(k,:)=allgenebins{gene_index,k}(end:-1:1);
end;

% plot original profiles without normalization
figure;
plotprofiles(ensemblid,tempprofile,1,10);

% apply normalization of time points (geommeanmedian, based on last time points)
normalized_tempprofile1=[];
for k=1:10,
  normalized_tempprofile1(k,:)=tempprofile(k,:)*normalizernumbers(k);
end;
plotprofiles(ensemblid,normalized_tempprofile,1,10);

% compute mean profile
meanprofile=mean(normalized_tempprofile1,1);

% apply normalization by mean profile
normalized_tempprofile2=[];
for k=1:10,
  normalized_tempprofile2(k,:)=normalized_tempprofile1(k,:)./meanprofile;
end;
plotprofiles(ensemblid,normalized_tempprofile2,1,10);

figure;
plotprofiles(ensemblid,meanprofile,1,1);

% compute center of mass, first without normalization with the mean profile
centerofmass1=[];
for k=1:10,
  centerofmass1(k)=sum(normalized_tempprofile1(k,:).*[1:size(normalized_tempprofile1,2)]*200)/...
    sum(normalized_tempprofile1(k,:));
end;
plotprofiles(ensemblid,normalized_tempprofile1,1,10,centerofmass1);

% compute center of mass, using normalization with the mean profile
centerofmass2=[];
for k=1:10,
centerofmass2(k)=sum(normalized_tempprofile2(k,:).*[1:size(normalized_tempprofile2,2)]*200)/...
sum(normalized_tempprofile2(k,:));
end;
plotprofiles(ensemblid,normalized_tempprofile2,1,10,centerofmass2);


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

compute_allcentersofmass;

normalized_masscenters1=masscenters1(:,2:6);
for m=1:size(masscenters1,1),
  m
  normalized_masscenters1(m,:)=(normalized_masscenters1(m,:)-mean(normalized_masscenters1(m,:)))/sqrt(var(normalized_masscenters1(m,:)));
end;
normalized_masscenters2=masscenters2(:,2:6);
  for m=1:size(masscenters2,1),
  m
  normalized_masscenters2(m,:)=(normalized_masscenters2(m,:)-mean(normalized_masscenters2(m,:)))/sqrt(var(normalized_masscenters2(m,:)));
end;


I=find(~isnan(masscenters2(:,1)));


nclusters=9;
initindices=ceil(length(I)*rand(nclusters,1));
options(1)=0;
options(2)=1e-5;
options(3)=1e-5;
options(14)=2000;
allclustindices=nan*ones(size(masscenters1,1),nclusters);
[clustmeans, clustoptions, clustindices, errlog]=kmeans(normalized_masscenters1(I(initindices),:), normalized_masscenters1(I,:), options);
allclustindices(I,:)=clustindices;
clustsizes=zeros(nclusters,1);
clustvariances=zeros(nclusters,5);
for k=1:nclusters,
  clustsizes(k)=sum(allclustindices(:,k)==1);

  Itemp=find(allclustindices(:,k)==1);
  clustvariances(k,:)=var(normalized_masscenters1(Itemp,:));
end;


[y,Isize]=sort(-clustsizes);
figure;
ymin=-2;
ymax=2;
for clustindex=1:9,
  h=subplot(3,3,clustindex);
  h=plot([5 10 20 40 80],clustmeans(Isize(clustindex),:),'b-.');
hold on;
  set(h,'linewidth',3);
h=plot([5 10 20 40 80],clustmeans(Isize(clustindex),:)-sqrt(clustvariances(Isize(clustindex),:)),'r--');
h=plot([5 10 20 40 80],clustmeans(Isize(clustindex),:)+sqrt(clustvariances(Isize(clustindex),:)),'r--');
  h=gca;
  set(h,'ylim',[ymin ymax]);
  title(sprintf('cluster %d, %d genes', clustindex, clustsizes(Isize(clustindex))));
end;

clustmeans1=clustmeans;
clustindices1=allclustindices;


nclusters=9;
initindices=ceil(length(I)*rand(nclusters,1));
options(1)=0;
options(2)=1e-5;
options(3)=1e-5;
options(14)=2000;
allclustindices=nan*ones(size(masscenters2,1),nclusters);
[clustmeans, clustoptions, clustindices, errlog]=kmeans(normalized_masscenters2(I(initindices),:), normalized_masscenters2(I,:), options);
allclustindices(I,:)=clustindices;
clustsizes=zeros(nclusters,1);
clustvariances=zeros(nclusters,5);
for k=1:nclusters,
clustsizes(k)=sum(allclustindices(:,k)==1);

Itemp=find(allclustindices(:,k)==1);
clustvariances(k,:)=var(normalized_masscenters2(Itemp,:));
end;


[y,Isize]=sort(-clustsizes);
figure;
ymin=-2;
ymax=2;
for clustindex=1:9,
h=subplot(3,3,clustindex);
h=plot([5 10 20 40 80],clustmeans(Isize(clustindex),:),'b-.');
hold on;
set(h,'linewidth',3);
h=plot([5 10 20 40 80],clustmeans(Isize(clustindex),:)-sqrt(clustvariances(Isize(clustindex),:)),'r--');
h=plot([5 10 20 40 80],clustmeans(Isize(clustindex),:)+sqrt(clustvariances(Isize(clustindex),:)),'r--');
h=gca;
set(h,'ylim',[ymin ymax]);
title(sprintf('cluster %d, %d genes', clustindex, clustsizes(Isize(clustindex))));
end;


normalizernumbers=[];
greb1=[];

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_last_0kb_to_4kb;
pol_summaryseries_normalized=normalizedbyfilteredcount_pol_summaryseries_last_0kb_to_4kb;
normalizernumbers(1,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(1,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_last_0kb_to_4kb;
pol_summaryseries_normalized=normalizedbytotalreadcount_pol_summaryseries_last_0kb_to_4kb;
normalizernumbers(2,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(2,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_last_0kb_to_4kb;
pol_summaryseries_normalized=normalizedbygeommeanmedian_pol_summaryseries_last_0kb_to_4kb;
normalizernumbers(3,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(3,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_last_0kb_to_4kb;
pol_summaryseries_normalized=normalizedbynoiselevel_pol_summaryseries_last_0kb_to_4kb;
normalizernumbers(4,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(4,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_first_0kb_to_4kb;
pol_summaryseries_normalized=normalizedbyfilteredcount_pol_summaryseries_first_0kb_to_4kb;
normalizernumbers(5,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(5,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_first_0kb_to_4kb;
pol_summaryseries_normalized=normalizedbytotalreadcount_pol_summaryseries_first_0kb_to_4kb;
normalizernumbers(6,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(6,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_first_0kb_to_4kb;
pol_summaryseries_normalized=normalizedbygeommeanmedian_pol_summaryseries_first_0kb_to_4kb;
normalizernumbers(7,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(7,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_first_0kb_to_4kb;
pol_summaryseries_normalized=normalizedbynoiselevel_pol_summaryseries_first_0kb_to_4kb;
normalizernumbers(8,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(8,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_first_1kb_to_5kb;
pol_summaryseries_normalized=normalizedbyfilteredcount_pol_summaryseries_first_1kb_to_5kb;
normalizernumbers(9,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(9,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_first_1kb_to_5kb;
pol_summaryseries_normalized=normalizedbytotalreadcount_pol_summaryseries_first_1kb_to_5kb;
normalizernumbers(10,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(10,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_first_1kb_to_5kb;
pol_summaryseries_normalized=normalizedbygeommeanmedian_pol_summaryseries_first_1kb_to_5kb;
normalizernumbers(11,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(11,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_first_1kb_to_5kb;
pol_summaryseries_normalized=normalizedbynoiselevel_pol_summaryseries_first_1kb_to_5kb;
normalizernumbers(12,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(12,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_meanovergene;
pol_summaryseries_normalized=normalizedbyfilteredcount_pol_summaryseries_meanovergene;
normalizernumbers(13,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(13,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_meanovergene;
pol_summaryseries_normalized=normalizedbytotalreadcount_pol_summaryseries_meanovergene;
normalizernumbers(14,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(14,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_meanovergene;
pol_summaryseries_normalized=normalizedbygeommeanmedian_pol_summaryseries_meanovergene;
normalizernumbers(15,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(15,:)=pol_summaryseries_normalized(6586,:);

pol_summaryseries_unnormalized=unnormalized_pol_summaryseries_meanovergene;
pol_summaryseries_normalized=normalizedbynoiselevel_pol_summaryseries_meanovergene;
normalizernumbers(16,:)=pol_summaryseries_normalized(2000,:)./pol_summaryseries_unnormalized(2000,:);
greb1(16,:)=pol_summaryseries_normalized(6586,:);

for k=1:16,normalizernumbers(k,:)=normalizernumbers(k,:)/mean(normalizernumbers(k,:));end;

for k=1:16,greb1(k,:)=greb1(k,:)/mean(greb1(k,:));end;



pol_summaryseries=normalizedbygeommeanmedian_pol_summaryseries_last_0kb_to_4kb;

pol_summaryseries_means=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
pol_summaryseries_means(k)=mean(pol_summaryseries(k,:));
end;
pol_summaryseries_variances=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
pol_summaryseries_variances(k)=var(pol_summaryseries(k,:));
end;
rna_summaryseries_means=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
rna_summaryseries_means(k)=mean(rna_summaryseries(k,:));
end;
rna_summaryseries_variances=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
rna_summaryseries_variances(k)=var(rna_summaryseries(k,:));
end;

% Isubstantialpol2=find((pol_summaryseries_means>=exp(9)) & ...
% (pol_summaryseries_variances>=exp(15)) & ...
% (rna_filteringresults==1));

Isubstantialpol2=find((pol_summaryseries_means>=exp(13.5)) & ...
(pol_summaryseries_variances>=exp(24)) & ...
(rna_filteringresults==1));


for k=1:size(pol_summaryseries,1),
pol_summaryseries(k,:)=(pol_summaryseries(k,:)-mean(pol_summaryseries(k,:)))/sqrt(var(pol_summaryseries(k,:)));
end;

for k=1:size(pol_summaryseries,1),
rna_summaryseries(k,:)=(rna_summaryseries(k,:)-mean(rna_summaryseries(k,:)))/sqrt(var(rna_summaryseries(k,:)));
end;

pol_summaryseries_backup=pol_summaryseries;
pol_summaryseries=[pol_summaryseries rna_summaryseries];

nclusters=9;
initindices=ceil(length(Isubstantialpol2)*rand(nclusters,1));
options(1)=0;
options(2)=1e-5;
options(3)=1e-5;
options(14)=2000;
allclustindices=nan*ones(size(pol_summaryseries,1),nclusters);
[clustmeans, clustoptions, clustindices, errlog]=kmeans(pol_summaryseries(Isubstantialpol2(initindices),:), pol_summaryseries(Isubstantialpol2,:), options);
allclustindices(Isubstantialpol2,:)=clustindices;
clustsizes=zeros(nclusters,1);
%clustvariances=zeros(nclusters,20);
clustvariances=zeros(nclusters,size(pol_summaryseries,2));
for k=1:nclusters,
clustsizes(k)=sum(allclustindices(:,k)==1);

Itemp=find(allclustindices(:,k)==1);
clustvariances(k,:)=var(pol_summaryseries(Itemp,:));
end;


[y,Isize]=sort(-clustsizes);
figure;
%ymin=-2;
%ymax=2;
for clustindex=1:9,
 % timescale=[0 5 10 20 40 80 160 320 640 1280];
  timescale=[1:10];
  h=subplot(3,3,clustindex);
  h=plot(timescale,clustmeans(Isize(clustindex),1:10),'b-.');
  hold on;
  set(h,'linewidth',3);
  h=plot(timescale,clustmeans(Isize(clustindex),1:10)-sqrt(clustvariances(Isize(clustindex),1:10)),'c--');
  h=plot(timescale,clustmeans(Isize(clustindex),1:10)+sqrt(clustvariances(Isize(clustindex),1:10)),'c--');
  if size(clustmeans,2)>10,
    h=plot(timescale,clustmeans(Isize(clustindex),11:20),'r-.');
    hold on;
    set(h,'linewidth',3);
    h=plot(timescale,clustmeans(Isize(clustindex),11:20)-sqrt(clustvariances(Isize(clustindex),11:20)),'m--');
    h=plot(timescale,clustmeans(Isize(clustindex),11:20)+sqrt(clustvariances(Isize(clustindex),11:20)),'m--');
  end;  
  h=gca;
  %set(h,'ylim',[ymin ymax]);
  title(sprintf('newest cluster %d, %d genes', clustindex, clustsizes(Isize(clustindex))));
end;







