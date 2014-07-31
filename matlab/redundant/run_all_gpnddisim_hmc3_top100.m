%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

%mybasedir_code='/share/work/jtpelto/tempsynergy/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
%mybasedir_code='~/synergy_data/tempcodebranch/';
%mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
mybasedir_code='~/mlprojects/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/share/work/jtpelto/synergy-data/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='~/synergy_data/';
mybasedir_data='~/projects/pol2rnaseq/';

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


load([mybasedir_data, 'data/series_for_matti_corrected_ver1.mat']);


good_genes = { 'ENSG00000152413',
'ENSG00000175901',
'ENSG00000186628',
'ENSG00000217791',
'ENSG00000178401',
'ENSG00000230807',
'ENSG00000060303',
'ENSG00000235036',
'ENSG00000230207',
'ENSG00000249822',
'ENSG00000242013',
'ENSG00000247017',
'ENSG00000198042',
'ENSG00000152804',
'ENSG00000035403',
'ENSG00000198744',
'ENSG00000165863',
'ENSG00000105976',
'ENSG00000173542',
'ENSG00000135940',
'ENSG00000113369',
'ENSG00000233287',
'ENSG00000125962',
'ENSG00000120798',
'ENSG00000225920',
'ENSG00000165185',
'ENSG00000243404',
'ENSG00000133048',
'ENSG00000251536',
'ENSG00000128872',
'ENSG00000196678',
'ENSG00000213297',
'ENSG00000028310',
'ENSG00000234841',
'ENSG00000088325',
'ENSG00000134533',
'ENSG00000081665',
'ENSG00000019102',
'ENSG00000231610',
'ENSG00000154040',
'ENSG00000204308',
'ENSG00000152253',
'ENSG00000213073',
'ENSG00000073282',
'ENSG00000107562',
'ENSG00000251948',
'ENSG00000182718',
'ENSG00000186815',
'ENSG00000249660',
'ENSG00000140987',
'ENSG00000175602',
'ENSG00000125735',
'ENSG00000160813',
'ENSG00000179021',
'ENSG00000162852',
'ENSG00000249921',
'ENSG00000122965',
'ENSG00000164983',
'ENSG00000180769',
'ENSG00000138468',
'ENSG00000164463',
'ENSG00000105379',
'ENSG00000142178',
'ENSG00000139514',
'ENSG00000183808',
'ENSG00000226745',
'ENSG00000168300',
'ENSG00000181444',
'ENSG00000133639',
'ENSG00000173402',
'ENSG00000115306',
'ENSG00000116001',
'ENSG00000130649',
'ENSG00000168702',
'ENSG00000137992',
'ENSG00000175895',
'ENSG00000054965',
'ENSG00000163362',
'ENSG00000134897',
'ENSG00000113312',
'ENSG00000119321',
'ENSG00000183161',
'ENSG00000249969',
'ENSG00000187735',
'ENSG00000142207',
'ENSG00000068784',
'ENSG00000196843',
'ENSG00000196659',
'ENSG00000125652',
'ENSG00000052749',
'ENSG00000122779',
'ENSG00000085231',
'ENSG00000198691',
'ENSG00000213867',
'ENSG00000177728',
'ENSG00000138735',
'ENSG00000168101',
'ENSG00000241279',
'ENSG00000180953',
'ENSG00000163655' };

indices = zeros(size(good_genes));
for k=1:length(good_genes),
  t = good_genes{k};
  indices(k) = find(str2double(t(10:end)) == bininfo(:, 5));
end





%---------------------------------------------------
% Load maximum likelihood fitting results with small and long delays
%---------------------------------------------------
cd(analysisdir)

%load allresults_shifted_polorrna4.mat
% allresults_ensemblids_pol2 allresults_geneindices_pol2 allresults_jointmodels_pol2 allresults_jointtransforminfos_pol2 allresults_loglikelihoods_pol2
% allresults_ensemblids_rna allresults_geneindices_rna allresults_jointmodels_rna allresults_jointtransforminfos_rna allresults_loglikelihoods_rna
%allresults_jointtransforminfos_pol2 = allresults_tinfos_pol2;
%allresults_jointtransforminfos_rna = allresults_tinfos_rna;


%load allresults_shifted_longerdelay5.mat
%allresults_ensemblids_joint = allresults_ensemblids; 
%allresults_geneindices_joint = allresults_geneindices;
%allresults_jointmodels_joint = allresults_jointmodels;
%allresults_jointtransforminfos_joint = allresults_jointtransforminfos;
%allresults_loglikelihoods_joint = allresults_loglikelihoods;




%---------------------------------------------------
% Compute estimated probability of joint fit being useful
%---------------------------------------------------

% probcomparison = nan*ones(max([length(allresults_jointtransforminfos_pol2) ...
%      length(allresults_jointtransforminfos_joint)]),1);
% for i=1:length(probcomparison),
%   if (i<=size(allresults_loglikelihoods_joint,1)) ...
%     && (i<=size(allresults_loglikelihoods_pol2,1)),
%     if (~isempty(allresults_jointmodels_joint{i})) && (~isempty(allresults_jointmodels_pol2{i})),
%       probcomparison(i)=allresults_loglikelihoods_joint(i,3) ...
%           -allresults_loglikelihoods_pol2(i,3) -allresults_loglikelihoods_rna(i,3);
%     end;
%   end;  
% end;


% Rank genes by probability of joint fit being useful
%[y,Ijointrank]=sort(-probcomparison);


%maxHMCgenes = 10;
nHMCiters = 10000;
%HMCsamples = cell(length(allresults_jointmodels_joint),1);

%istart=1;
%iend=maxHMCgenes;
for i=myI,
  % Start sampling from the maximum likelihood joint-model fit
  %model = allresults_jointmodels_joint{Ijointrank(i)};

  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  %model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  %model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  %model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;
  gene_index = indices(i);
  gene_name = good_genes{i};
  
  randn('seed',bininfo(gene_index,5));
  rand('seed',bininfo(gene_index,5)+1234567);
  
  dataVals1=pol_summaryseries(gene_index,:)';
  dataVals2=rna_summaryseries(gene_index,:)';
  timevector=[0 5 10 20 40 80 160 320 640 1280]' + timeshift;

  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;
    
  lengthscale=2;

  timevector = {timevector,timevector};
  dataVals = {tempvals1, tempvals2};
  initializationtype=1;

  [m,temptransforminfo]=createNdSimDisim_celltimes(timevector, ...
						   dataVals,lengthscale,initializationtype,[],[],1);

  % Apply HMC sampling
  HMCsamples = gpnddisimSampleHMC(m, 1, nHMCiters);
  HMCsamples = HMCsamples(10:10:end, :);
  save(sprintf('hmc_results/%s_samples_%s.mat', gene_name, id), ...
       'gene_name', 'gene_index', 'm', 'HMCsamples');
end;

