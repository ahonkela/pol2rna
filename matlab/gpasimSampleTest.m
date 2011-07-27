basedir='~/mlprojects/';
pol2dir='~/synergy_data/PolII/Mapping_results/'

% for model-level computations
addpath([basedir 'gpsim/matlab/jaakko_testversion/']);

% for kernel-level computations
addpath([basedir 'kern/matlab/jaakko_testversion/']);

% for optimiDefaultConstraint.m
addpath([basedir 'optimi/matlab/']);

% for lnDiffErfs.m
addpath([basedir 'ndlutil/matlab/']);

% for addPrior.m
addpath([basedir 'prior/matlab/']);

% for dist2.m
addpath([basedir 'matlab/netlab/NETLAB3p3/']);

% for modelTieParam.m
addpath([basedir 'mltools/matlab/']);

cd(pol2dir)
load pol2_for_matti_ver2.mat
pol_summaryseries=normalizedbygeommeanmedian_pol_summaryseries_last_0kb_to_4kb;


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
options=struct();
options.includeNoise=1;
options.addPriors=1;
options.asynchronyType=1;
options.optimiser='MH'; % Metropolis-Hastings sampler
annotation=[];

model = gpasimCreate(numGenes, times, dataVals, dataVars, options, annotation );


dataVals{1}=randn(10,1);
model = gpasimCreate(0, {times{1}}, {dataVals{1}}, {dataVars{1}}, options, annotation );

dim = size(dataVals{1}, 1);
tempK=eye(dim)*var(dataVals{1});
tempm=dataVals{1}-mean(dataVals{1});
ll = -dim*log(2*pi) - log(det(tempK)) - tempm'*inv(tempK)*tempm;
ll = ll*0.5



display=1;
iters=10000;
proposalvariance=0.1;
proposalvariance_polmean=0.000002*var(dataVals{1});
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
