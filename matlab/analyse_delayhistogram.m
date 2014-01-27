%mybasedir_code='/share/work/jtpelto/tempsynergy/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
mybasedir_code='~/synergy_data/tempcodebranch/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='~/mlprojects/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
mybasedir_data='~/synergy_data/';

mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';


% pol2 and H3K4me3 data
pol2dir=[mybasedir_data 'PolII/Mapping_results/'];
h3k4me3dir=[mybasedir_data 'H3K4me3/Mapping_results/'];
analysisdir=[mybasedir_analyses 'analyses/'];

% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/jaakko_testversion/']
% for model-level computations
path2=[mybasedir_code 'gpsim/matlab/jaakko_testversion/'];
% for optimiDefaultConstraint.m
path3=[mybasedir_code 'optimi/matlab/jaakko_testversion/'];
% for lnDiffErfs.m
path4=[mybasedir_code 'ndlutil/matlab/'];
% for addPrior.m
path5=[mybasedir_code 'prior/matlab/'];
% for dist2.m
path6=[mybasedir_code 'matlab/netlab/NETLAB3p3/'];
% for modelTieParam.m
path7=[mybasedir_code 'mltools/matlab/'];
% for various experiment things
path8=[mybasedir_code 'pol2rnaseq/matlab/'];

addpath(path1,path2,path3,path4,path5,path6,path7,path8)



%---------------------------------------------------
% Load gene data
%---------------------------------------------------
cd(h3k4me3dir)
load series_for_matti_corrected_ver1.mat  % provides pol_summaryseries,rna_summaryseries

cd(analysisdir)
load models_Pol2.mat % provides interestinggenes
templls=lls_Pol2(:,1)-lls_Pol2(:,3);
[y,I]=sort(-templls);
interestinggenes=genes_Pol2;
interestinggenes=interestinggenes(I);


%---------------------------------------------------
% Find significant genes according to GPREGE ranking
%---------------------------------------------------
cd(analysisdir)

rnarankingfile=read_stringfile('gpregeranking_rna.txt', [32 9 10 13]);
rnarankingscores=ones(length(interestinggenes),1)*nan;
for k=1:length(rnarankingfile)-1,
  rnarankingscores(k)=str2double(rnarankingfile{k+1}{2});
end;
rnarankingprobs=rnarankingscores./(rnarankingscores+1);

pol2rankingfile=read_stringfile('gpregeranking_pol2.txt', [32 9 10 13]);
pol2rankingscores=ones(length(interestinggenes),1)*nan;
for k=1:length(rnarankingfile)-1,
  pol2rankingscores(k)=str2double(pol2rankingfile{k+1}{2});
end;
pol2rankingprobs=pol2rankingscores./(pol2rankingscores+1);

pol2vars=var(pol_summaryseries(interestinggenes,:)')';

significantgenes=find((rnarankingprobs>0.70) & (pol2rankingprobs>0.70) & (pol2vars>0));
significantgenes=interestinggenes(significantgenes);


%---------------------------------------------------
% Compute the delays between peaks
%---------------------------------------------------
times=[0 5 10 20 40 80 160 320 640 1280];
delays=ones(length(significantgenes),1)*nan;
for k=1:length(significantgenes),
  gene_index=significantgenes(k);
  Ipol2=find(pol_summaryseries(gene_index,:)==max(pol_summaryseries(gene_index,:)));
  Irna=find(rna_summaryseries(gene_index,:)==max(rna_summaryseries(gene_index,:)));
  delays(k)=times(Irna)-times(Ipol2);
end;


%---------------------------------------------------
% Plot the series
%---------------------------------------------------
k=1;while k<length(significantgenes),
  gene_index=significantgenes(k);
  tempseries=pol_summaryseries(gene_index,:);
  tempseries=(tempseries-min(tempseries))./(max(tempseries)-min(tempseries));
  
  clf;
  plot(sqrt([0 5 10 20 40 80 160 320 640 1280]), ...
       tempseries,'r-');
  title(sprintf('%d',k));
  hold on;
  tempseries=rna_summaryseries(gene_index,:);
  tempseries=(tempseries-min(tempseries))./(max(tempseries)-min(tempseries));
  plot(sqrt([0 5 10 20 40 80 160 320 640 1280]), ...
       tempseries,'g-');  
  
  x=kbhit();
  if x=='w',k=k-1;
  else k=k+1;end;
end;

