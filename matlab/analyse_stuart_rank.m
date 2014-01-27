% Analyse how high is a list of genes provided by Stuart Wilson, 
% whose website is at http://www.shef.ac.uk/mbb/staff/wilson
% The list of genes is currently confidential.
% Idea: for each of Stuart's genes, plot a position in a 2D scatterplot,
%   where the coordinates are: the rank of the gene in our "usefulness 
%   of joint GP modeling" ranking, and the Pol2-RNA delay estimated by 
%   the joint GP fit.
% If our estimated delay for Stuart's genes is typically high, and 
% conversely if many of our genes appear in one of Stuart's lists,
% this can indicate that Stuart's explanation for Pol2-RNA delay could
% apply to our data set.
%
% Jaakko Peltonen, Nov 2, 2011.

%mybasedir_code='/media/JPELTONEN4/mlprojects/';
mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='~/mlprojects/';
%mybasedir_code='/share/work/jtpelto/tempsynergy/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='~/synergy_data/';
mybasedir_data='/share/work/jtpelto/synergy-data/';

mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';

% pol2 and H3K4me3 data
pol2dir=[mybasedir_data 'PolII/Mapping_results/'];
h3k4me3dir=[mybasedir_data 'H3K4me3/Mapping_results/'];
analysisdir=[mybasedir_analyses 'analyses/'];

% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/']
% for model-level computations
path2=[mybasedir_code 'gpsim/matlab/'];
% for optimiDefaultConstraint.m
path3=[mybasedir_code 'optimi/matlab/'];
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
% Load fitting results with small and long delays
%---------------------------------------------------
cd(analysisdir)

load allresults_shifted_polorrna4.mat
% allresults_ensemblids_pol2 allresults_geneindices_pol2 allresults_jointmodels_pol2 allresults_jointtransforminfos_pol2 allresults_loglikelihoods_pol2
% allresults_ensemblids_rna allresults_geneindices_rna allresults_jointmodels_rna allresults_jointtransforminfos_rna allresults_loglikelihoods_rna
allresults_jointtransforminfos_pol2 = allresults_tinfos_pol2;
allresults_jointtransforminfos_rna = allresults_tinfos_rna;

load allresults_shifted_longerdelay5.mat
allresults_ensemblids_joint = allresults_ensemblids; 
allresults_geneindices_joint = allresults_geneindices;
allresults_jointmodels_joint = allresults_jointmodels;
allresults_jointtransforminfos_joint = allresults_jointtransforminfos;
allresults_loglikelihoods_joint = allresults_loglikelihoods;

%---------------------------------------------------
% Compute estimated probability of gene having long delay
%---------------------------------------------------

probcomparison = nan*ones(max([length(allresults_jointtransforminfos_pol2) ...
     length(allresults_jointtransforminfos_joint)]),1);
for i=1:length(probcomparison),
  if (i<=size(allresults_loglikelihoods_joint,1)) ...
    && (i<=size(allresults_loglikelihoods_pol2,1)),
    if (~isempty(allresults_jointmodels_joint{i})) && (~isempty(allresults_jointmodels_pol2{i})),
      probcomparison(i)=allresults_loglikelihoods_joint(i,3) ...
          -allresults_loglikelihoods_pol2(i,3) -allresults_loglikelihoods_rna(i,3);
    end;
  end;  
end;

[y,Iprobcomparison]=sort(-probcomparison);




% Load translation table from ENSEMBL transcript ids to ENSEMBL gene ids
% (this is needed because stuart uses transcript ids, we use ENSEMBL gene ids)
stuartdir='/share/work/jtpelto/synergy-data/analyses/stuartgenes/Stuart/';
if 0,
  filename='ensembl_transcript_to_gene.txt';
  translationtable=read_stringfile([stuartdir filename],[' ' 10 13 0],'');
  % first column of translationmatrix is gene id, second column is transcript id
  translationmatrix=nan*ones(length(translationtable),2);
  for k=1:length(translationtable),
    if mod(k,100)==0, fprintf(1,'%d\n',k); end;
    if strcmp(translationtable{k}{1}(1:5),'"ENSG') == 1,
      translationmatrix(k,1)=str2double(translationtable{k}{1}(6:end-1));
    end;
    if strcmp(translationtable{k}{2}(1:5),'"ENST') == 1,
      translationmatrix(k,2)=str2double(translationtable{k}{2}(6:end-1));
    end;
  end;
  save([stuartdir 'ensembl_transcript_to_gene.mat'],'translationmatrix','-mat');
end;
load([stuartdir 'ensembl_transcript_to_gene.mat'],'translationmatrix','-mat');


% Load Stuart's gene list
filenames=cell(11,1); idcolumns=zeros(11,1);
filenames{1}='CYTO_DDX39_DOWN_FC_1.5.csv'; idcolumns(1)=4;
filenames{2}='CYTO_SRAG_REF_DOWN_FC_1.5.csv'; idcolumns(2)=4;
filenames{3}='CYTO_UAP56_DDX39_DOWN_FC_1.5.csv'; idcolumns(3)=4;
filenames{4}='CYTO_REF_DOWN_FC_1.5.csv'; idcolumns(4)=4;
filenames{5}='CYTO_SSRP1_DOWN_FC_1.5.csv'; idcolumns(5)=4;
filenames{6}='CYTO_UAP56_DOWN_FC_1.5.csv'; idcolumns(6)=4;
filenames{7}='CYTO_SPT16_DOWN_FC_1.5.csv'; idcolumns(7)=4;
filenames{8}='CYTO_THOC5_DOWN_FC_1.5.csv'; idcolumns(8)=4;
filenames{9}='CYTO_UIF_DOWN_FC_1.5.csv'; idcolumns(9)=4;
filenames{10}='CYTO_SRAG_DOWN_FC_1.5.csv'; idcolumns(10)=4;
filenames{11}='CYTO_THOC5_REF_DOWN_FC_1.5.csv'; idcolumns(11)=5;



stuartranks=nan*ones(length(allresults_jointtransforminfos_joint),11);

for stuartfileindex=1:11,
  filename = filenames{stuartfileindex};
  idcolumn = idcolumns(stuartfileindex);
  fprintf(1,'Processing Stuart file %s\n',filenames{stuartfileindex});
  stuartfile=read_stringfile([stuartdir filename],'','#');

  % extract transcript ids from Stuart's file
  stuart_geneids=zeros(length(stuartfile)-1,1);
  for k=1:length(stuart_geneids),
    stuarttext=stuartfile{k+1}{idcolumn}(6:end-1);
    if length(stuarttext)>0,
      stuart_geneids(k) = str2double(stuarttext);
    else
      stuart_geneids(k) = nan;
    end;
  end;
  % translate transcript ids to gene ids
  stuart_geneids2=stuart_geneids;
  for k=1:length(stuart_geneids2),
    if mod(k,100)==0, fprintf(1,'%d\n',k); end;
    %fprintf(1,'%d\n',k);
    I=find(translationmatrix(:,2)==stuart_geneids(k),1);
    if length(I)==1,
      stuart_geneids2(k)=translationmatrix(I(1),1);
    elseif length(I)>1,
      fprintf(1,'Error, found multiple transcripts with same id at k=%d\n',k);
    end;
  end;

  % Assign ranks from stuart's list to our genes
  for k=1:length(allresults_jointtransforminfos_pol2),
    I=find(allresults_ensemblids_joint(k)==stuart_geneids2);
    if length(I)==1,
      stuartranks(k,stuartfileindex)=I(1);
    elseif length(I)>1,
      stuartranks(k,stuartfileindex)=I(1);
      fprintf(1,'Warning, found multiple stuart-genes with same id at k=%d\n',k);
    end;     
  end;

end;



genedelays=nan*ones(length(probcomparison),1);
for i=1:length(probcomparison),
  if (~isempty(allresults_jointmodels_joint{i}))
    genedelays(i) = allresults_jointmodels_joint{i}.delay;
  end;
end;




% Proportion of our genes that have reasonably large delays
delaycutoff=30;
sum(genedelays>=delaycutoff)/sum(~isnan(genedelays))

% Number of our genes on Stuart's various lists
a1=sum((~isnan(stuartranks)),1)

% Number of our genes that have reasonably large delays on Stuart's various lists
a2=sum(repmat(genedelays>=delaycutoff,[1 size(stuartranks,2)]).*(~isnan(stuartranks)),1)

% Proportion of the above: does being on Stuart's list generally mean larger delay
a2./a1



% Genes that are well modeled by our fits
wellmodeled = zeros(size(genedelays));
wellmodeled(Iprobcomparison(1:1000))=1;

% Proportion of well-fit genes that have reasonably large delays
delaycutoff=30;
sum((genedelays>=delaycutoff)&(wellmodeled==1))/sum(wellmodeled==1)

% Number of our genes on Stuart's various lists
a1=sum((~isnan(stuartranks))&repmat((wellmodeled==1),[1 size(stuartranks,2)]),1)

% Number of our genes that have reasonably large delays on Stuart's various lists
a2=sum(repmat(genedelays>=delaycutoff,[1 size(stuartranks,2)]).*repmat((wellmodeled==1),[1 size(stuartranks,2)]).*(~isnan(stuartranks)),1)

% Proportion of the above: does being on Stuart's list generally mean larger delay
a2./a1
