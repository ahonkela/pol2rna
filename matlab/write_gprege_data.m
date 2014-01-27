
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



cd(h3k4me3dir)
load series_for_matti_corrected_ver1.mat

cd(analysisdir)
load models_Pol2.mat
templls=lls_Pol2(:,1)-lls_Pol2(:,3);
[y,I]=sort(-templls);
interestinggenes=genes_Pol2;
interestinggenes=interestinggenes(I);



n_interesting_genes=length(interestinggenes);
n_interesting_genes
startpercent=0;
kstart=floor(startpercent*n_interesting_genes/100)+1;
kend=n_interesting_genes;


fpol2=fopen('series_for_matti_gprege_pol2.txt','w');
frna=fopen('series_for_matti_gprege_rna.txt','w');


timevector=[0 5 10 20 40 80 160 320 640 1280]';
for k=1:length(timevector), 
  fprintf(fpol2,'\"%d\"', timevector(k));
  if k==length(timevector), fprintf(fpol2,'\n');else fprintf(fpol2,' ');end;
  fprintf(frna,'\"%d\"', timevector(k));
  if k==length(timevector), fprintf(frna,'\n');else fprintf(frna,' ');end;
end;
for k=1:length(interestinggenes),
  gene_index=interestinggenes(k);
  gene_name=sprintf('ENSG%d',bininfo(gene_index,5));
  fprintf(fpol2,'\"%s\" ', gene_name);
  fprintf(frna,'\"%s\" ', gene_name);
  
  dataVals1=pol_summaryseries(gene_index,:)';
  dataVals2=rna_summaryseries(gene_index,:)';
  
  for l=1:length(timevector),
    fprintf(fpol2,'%f', dataVals1(l));
    if l==length(timevector), fprintf(fpol2,'\n');else fprintf(fpol2,' ');end;
    fprintf(frna,'%f', dataVals2(l));
    if l==length(timevector), fprintf(frna,'\n');else fprintf(frna,' ');end;
  end;
end;
fclose(fpol2);
fclose(frna);


