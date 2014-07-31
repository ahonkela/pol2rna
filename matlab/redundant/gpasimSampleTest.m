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
% for various experiment things
path8=[mybasedir 'pol2rnaseq/matlab/'];

addpath(path1,path2,path3,path4,path5,path6,path7,path8)


cd(pol2dir)
load pol2_for_matti_ver3.mat
pol_summaryseries=normalizedbygeommeanmedian_pol_summaryseries_last20percent;


gene_index=find(bininfo(:,5)==196208);  % GREB1, ENSEMBL-id 196208;


numGenes=1;

times=cell(1,2);

% POL2 observation times
times{1}=[0 5 10 20 40 80 160 320 640 1280]';
%times{1}=[0 1 2 3 4 5 6 7 8 9]';


% mRNA observation times
times{2}=[0 5 10 20 40 80 160 320 640 1280]';
%times{2}=[0 1 2 3 4 5 6 7 8 9]';

dataVals=cell(1,2);
dataVals{1}=pol_summaryseries(gene_index,:)';
dataVals{2}=rna_summaryseries(gene_index,:)';
dataVars{1}=nan*pol_summaryseries(gene_index,:)';
dataVars{2}=nan*rna_summaryseries(gene_index,:)';


% create simple test data
tstep=0.01;
testt=[0:tstep:9];
testvals1=1.45*sin(testt);
testB=0;
testD=0.1;
testS=sqrt(10)*31.69;
testvals2=testB/testD+testS*exp(-testD*testt).*cumsum(exp(testD*testt).*testvals1*tstep)*tstep;
testobst=testt(1:100:end);
testobs1=testvals1(1:100:end);
testobs2=testvals2(1:100:end);
figure;
plot(testt,testvals1,'b--');hold on;plot(testt,testvals2,'r--');
plot(testobst,testobs1,'b*');hold on;plot(testobst,testobs2,'r*');

model=createSimDisim(testobst',testobs1',testobs2',1);
modelb=gpdisimOptimise(model);
figure; plotpredictions(modelb,[0:0.1:9]',0,1);
[pars,nams]=gpdisimExtractParam(modelb)
pars(4)=log(testD);
pars(5)=log(testS);
modelc=gpdisimExpandParam(modelb,pars);
figure; plotpredictions(modelc,[0:0.1:9]',0,1);

predtimes=[0:9]';
K_new=kernCompute(modelc.kern, predtimes);
K_new_old=kernCompute(modelc.kern, predtimes, model.t);
K_old=modelc.K;
K_old_new=K_new_old';
figure;imagesc(K_new_old*inv(K_old));


options=struct();
options.includeNoise=1;
options.addPriors=1;
options.asynchronyType=1;
%options.optimiser='MH'; % Metropolis-Hastings sampler
annotation=[];


% For comparison purposes, create a model based on the disim code
dataValsMatrix=[dataVals{1} dataVals{2}];
gpsimoptions=options;
gpsimoptions.optimiser='conjgrad';
gpsimmodel = gpasimTempCreate(1, 1, times{1}, dataValsMatrix, 0*dataValsMatrix, gpsimoptions, annotation )
tempfix.index=2;
tempfix.value=0;
gpsimmodel.fix=tempfix;
[gpsimparams,gpsimnames]=gpasimTempExtractParam(gpsimmodel)
gpsimmodel=gpasimTempExpandParam(gpsimmodel,gpsimparams);
gpasimTempGradient(gpasimTempExtractParam(gpsimmodel),gpsimmodel)
gpsimmodel = gpasimTempOptimise(gpsimmodel);


% For comparison purposes, create a SIM-DISIM model with decay 0
% for the SIM part
gpsim3options=options;
gpsim3options.optimiser='conjgrad';
includerna=0;
pol2issim=0;
pol2isrbf=1;
if (includerna==1) & (pol2issim==1),
  numgenes=1;
  dataValsMatrix=[dataVals{1} dataVals{2}];
  gpsimoptions.fix=struct();
  % decay of f to p, here fixed to 0
  gpsim3options.fix.index(1)=1;    
  gpsim3options.fix.value(1)=-20;
  % rbf variance, unidentifiable, here fixed to 1
  gpsim3options.fix.index(2)=6;     
  gpsim3options.fix.value(2)=0;  
  dataValsMatrix(:,1)=dataValsMatrix(:,1)-dataValsMatrix(1,1);
  gpsim3model = gpasimTemp3Create(numgenes, 1, times{1}, dataValsMatrix, ...
                                  0*dataValsMatrix, gpsim3options, annotation )
  [pars,nams]=gpdisimExtractParam(gpsim3model)
  % initialization of inverse squared width
  pars(2)=log(1/(15^2));
  % initialization of POL2 effect variance
  pars(3)=10;%log(0.05*var(dataVals{1}));
  % initialization of RNA decay
  pars(4)=-2;
  % initialization of POL2 noise variance
  pars(7)=log(0.001*var(dataVals{1}));
  % initialization of RNA noise variance
  pars(8)=log(0.02*var(dataVals{2}));
  % initialization of RNA basal rate
  pars(9)=1;  
  gpsim3model=gpdisimExpandParam(gpsim3model,pars);
  gpsim3modelb = gpdisimOptimise(gpsim3model);
  
elseif pol2issim==1,
  numgenes=0;
  dataValsMatrix=[dataVals{1}];
  gpsimoptions.fix=struct();
  % decay of f to p, here fixed to -1e6
  gpsim3options.fix.index(1)=1;    
  gpsim3options.fix.value(1)=-20;

  dataValsMatrix(:,1)=dataValsMatrix(:,1)-dataValsMatrix(1,1);
  gpsim3model = gpasimTemp3Create(numgenes, 1, times{1}, dataValsMatrix, ...
                                  0*dataValsMatrix, gpsim3options, annotation )
  [pars,nams]=gpdisimExtractParam(gpsim3model)
  % initialization of inverse squared width
  pars(2)=log(1/(15^2));
  % initialization of POL2 effect variance
  pars(3)=20;
  % initialization of POL2 noise variance
  pars(4)=log(0.05*var(dataVals{1}));
elseif pol2isrbf==1,
  gpsim3options.fix=[];
  numgenes=0;
  dataValsMatrix=[dataVals{1}];
  dataValsMatrix(:,1)=dataValsMatrix(:,1)-dataValsMatrix(1,1);
  gpsim3model = gpasimTemp3Create(numgenes, 1, times{1}, dataValsMatrix, ...
                                  0*dataValsMatrix, gpsim3options, annotation )
  [pars,nams]=gpdisimExtractParam(gpsim3model)
  % initialization of inverse squared width
  pars(1)=log(1/(15^2));
  % initialization of POL2 effect variance
  pars(2)=20;
  % initialization of POL2 noise variance
  pars(3)=log(0.05*var(dataVals{1}));
end;
gpsim3model=gpdisimExpandParam(gpsim3model,pars);

gpsim3modelb = gpdisimOptimise(gpsim3model);

%tempfix.index=2;
%tempfix.value=0;
%gpsim2model.fix=tempfix;
[gpsim3params,gpsim3names]=gpdisimExtractParam(gpsim3model)
gpsim2model=gpasimTemp2ExpandParam(gpsim2model,gpsim2params);
gpasim2TempGradient(gpasimTemp2ExtractParam(gpsim2model),gpsim2model)
gpsim2model = gpasimTemp2Optimise(gpsim2model);



model = gpasimCreate(numGenes, times, dataVals, dataVars, options, annotation );
tempparams=gpasimExtractParam(model);
model = gpasimExpandParam(model,tempparams);
model = gpasimInitialize(model, times, dataVals);
% add random jitter to handcrafted initialization
tempparams=gpasimExtractParam(model);
for k=1:length(tempparams),tempparams(k)=tempparams(k)+5*randn;end;
model = gpasimExpandParam(model,tempparams);
model = gpasimOptimise(model);

predicttimes=[0:1:1280];
pol2predictions=nan*predicttimes;
pol2vars=nan*predicttimes;
rnapredictions=nan*predicttimes;
rnavars=nan*predicttimes;
for k=1:length(predicttimes),
  k
  [predmeans,predcov]=gpasimTemp3Predict(gpsim3model,predicttimes(k),predicttimes(k));
  pol2predictions(k)=predmeans(1);
  rnapredictions(k)=predmeans(2);
  pol2vars(k)=predcov(1,1);
  rnavars(k)=predcov(2,2);
end;

figure; 
subplot(2,1,1);
h=plot(times{1},dataVals{1},'k-o');
hold on; plot(predicttimes,pol2predictions,'b-');
hold on; plot(predicttimes,pol2predictions+pol2vars.^0.5,'r-');
hold on; plot(predicttimes,pol2predictions-pol2vars.^0.5,'r-');
axis([min(times{1}) max(times{1}) min(dataVals{1})-sqrt(var(dataVals{1})) max(dataVals{1})+sqrt(var(dataVals{1}))]);

%figure;
subplot(2,1,2);
plot(times{2},dataVals{2},'k-o');
hold on; plot(predicttimes,rnapredictions,'g-');
hold on; plot(predicttimes,rnapredictions+pol2vars.^0.5,'c-');
hold on; plot(predicttimes,rnapredictions-pol2vars.^0.5,'c-');
axis([min(times{2}) max(times{2}) min(dataVals{2})-sqrt(var(dataVals{2})) max(dataVals{2})+sqrt(var(dataVals{2}))]);


%--------------
% Independent ARBF models for both POL2 and RNA
%--------------
model = gpasimCreate(0, {times{1}}, {dataVals{1}}, {dataVars{1}}, options, annotation );
tempparams=gpasimExtractParam(model);
model = gpasimExpandParam(model,tempparams);
model = gpasimInitialize(model, times, {dataVals{1}});
% add random jitter to handcrafted initialization
tempparams=gpasimExtractParam(model);
for k=1:length(tempparams),tempparams(k)=tempparams(k)+5*randn;end;
model = gpasimExpandParam(model,tempparams);
model = gpasimOptimise(model);
model1=model;

model = gpasimCreate(0, {times{2}}, {dataVals{2}}, {dataVars{2}}, options, annotation );
tempparams=gpasimExtractParam(model);
model = gpasimExpandParam(model,tempparams);
model = gpasimInitialize(model, times, {dataVals{2}});
% add random jitter to handcrafted initialization
tempparams=gpasimExtractParam(model);
for k=1:length(tempparams),tempparams(k)=tempparams(k)+5*randn;end;
model = gpasimExpandParam(model,tempparams);
model = gpasimOptimise(model);
model2=model;



%dataVals{1}=randn(10,1);
model = gpasimCreate(0, {times{1}}, {dataVals{1}}, {dataVars{1}}, options, annotation );

% Naive log-likelihood of POL2
dim = size(dataVals{1}, 1);
tempK=eye(dim)*var(dataVals{1});
tempm=dataVals{1}-mean(dataVals{1});
ll1 = -dim*log(2*pi) - log(det(tempK)) - tempm'*inv(tempK)*tempm;
ll1 = ll1*0.5

% Naive log-likelihood of RNA
dim = size(dataVals{2}, 1);
tempK=eye(dim)*var(dataVals{2});
tempm=dataVals{2}-mean(dataVals{2});
ll2 = -dim*log(2*pi) - log(det(tempK)) - tempm'*inv(tempK)*tempm;
ll2 = ll2*0.5

% Naive log-likelihood of independently modeled POL2 and RNA
ll=ll1+ll2


display=1;
iters=1000;
linesearchiters=20;
proposalvariance=0.1;
proposalvariance_polmean=0.000002*var(dataVals{1});
samples = gpasimOptimizeGradient(model, dataVals, display, iters, linesearchiters,proposalvariance, proposalvariance_polmean);


samples = gpasimSampleMH(model, display, iters, proposalvariance, proposalvariance_polmean);


D_i=exp(randn);
D_j=exp(randn);
delta_i=0;
delta_j=0;
inverseWidth=exp(randn);
sigma=sqrt(2/inverseWidth);
asyncparam=exp(randn);

asyncparam=0.62
t1=exp(randn); t1=640;
t2=exp(randn); t2=320;



val1=0;
summult=0;
for a1=-0:1:1280,
  a1
  val1
  for a2=-0:1:1280,
    tempval=simComputeH(a1, a2, D_i, D_j, delta_i, delta_j, sigma);
    tempmult1=exp(-(a1-t1)^2/(2*asyncparam*t1))/sqrt(2*pi*asyncparam*t1);
    tempmult2=exp(-(a2-t2)^2/(2*asyncparam*t2))/sqrt(2*pi*asyncparam*t2);
    val1=val1+tempval*tempmult1*tempmult2;
    summult=summult+tempmult1*tempmult2;
  end;    
end;
val1=val1/summult



arbfVariance=1;
[s1a,hpart1a,s2a,hpart2a]=asimComputeH_valueonly(...
t2, t1, ...
    D_j, D_i, ...
    delta_j, delta_i, ...
    inverseWidth, arbfVariance, ...
    1, asyncparam, asyncparam)
val1b=s1a*exp(hpart1a)-s2a*exp(hpart2a)

val1c=asimComputeH_valueonly(...
t1, t2, ...
    D_i, D_j, ...
    delta_i, delta_j, ...
    inverseWidth, arbfVariance, ...
    1, asyncparam, asyncparam)

tempkern=model.kern.comp{1}.comp{2}; %ASIM kernel
tempkern.decay=D_i;
tempkern.delay=0;
tempkern.inverseWidth=inverseWidth;
tempkern.arbfVariance=arbfVariance;
tempkern.asynchronyParam1=asyncparam;

asimKernCompute(tempkern, t1, t2)


tempval=simComputeH(t1, t2, D_i, D_j, delta_i, delta_j, sigma)
val1b=asimComputeH_valueonly(t2,t1,D_j,D_i,delta_j,delta_i,inverseWidth,arbfVariance,1,asyncparam,asyncparam)



tempkern=model.kern.comp{1}.comp{2}; %ASIM kernel
D_i=
D_j=D_i;
inverseWidth=0.020114;
asyncparam=0.623676;
arbfVariance=1;
tempkern.decay=D_i;
tempkern.delay=0;
tempkern.inverseWidth=inverseWidth
tempkern.arbfVariance=arbfVariance;
tempkern.asynchronyParam1=asyncparam;
t1=640;t2=320;
asimKernCompute(tempkern, t1, t2)

%tempkern.asynchronyParam1=0
asimKernCompute(tempkern, times{2}, times{2})


D=3.800971;
inverseWidth=0.020114;
sensitivity=1;
arbfVariance=1;
asyncparam=0.623676;
t=640;

S=sensitivity;
L=sqrt(2/inverseWidth);
A=arbfVariance;
h=asyncparam*t;

[erfd1,logerfd1,signerfd1]=...
    myerfdiff((t+D*L*L/2+D*h)/sqrt(L*L+2*h),(D*L*L/2+D*2*h)/sqrt(L*L+4*h))
[erfd2,logerfd2,signerfd2]=...
    myerfdiff(D*L/2,(-t+D*L*L/2+D*h)/sqrt(L*L+2*h))

hpart1=exp(D*D*L*L/4-log(2*D)+D*D*h)*myerfdiff((t+D*L*L/2+D*h)/sqrt(L*L+2*h),(D*L*L/2+D*2*h)/sqrt(L*L+4*h))...
			     -exp(D*D*L*L/4-log(2*D)-2*D*t+D*D*h)*myerfdiff(D*L/2,(-t+D*L*L/2+D*h)/sqrt(L*L+2*h))
hpart1b=signerfd1*exp(D*D*L*L/4-log(2*D)+D*D*h + logerfd1) ...
	-signerfd2*exp(D*D*L*L/4-log(2*D)-2*D*t+D*D*h + logerfd2)

k=hpart1*2*S*S*A*L*sqrt(pi)/2


sigma=sqrt(2/inverseWidth);
simComputeH(t1, t2, D_i, D_j, delta_i, delta_j, sigma)




tempsimkern=kernCreate({t}, 'sim')
tempsimkern.variance=sensitivity*sensitivity;
tempsimkern.inverseWidth=inverseWidth;
tempsimkern.decay=D;
tempsimkern.delay=0;


val1=0;
summult=0;
for a1=-0:1:1280,
  a1
  val1
  for a2=-0:1:1280,
    tempval=simKernCompute(tempsimkern,t)
    %tempval=simComputeH(a1, a2, D_i, D_j, delta_i, delta_j, sigma);
    tempmult1=exp(-(a1-t1)^2/(2*asyncparam*t1))/sqrt(2*pi*asyncparam*t1);
    tempmult2=exp(-(a2-t2)^2/(2*asyncparam*t2))/sqrt(2*pi*asyncparam*t2);
    val1=val1+tempval*tempmult1*tempmult2;
    summult=summult+tempmult1*tempmult2;
  end;    
end;
val1=val1/summult



% h(1,1) = Inf = (1)*exp((6123.819111)+(-6130.784094))-(-1)*exp((1258.576231)+(-362.653035))
