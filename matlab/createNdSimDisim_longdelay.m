function [gpsim3model,transforminfo]=createNdSimDisim_longdelay(timevector,dataVals1,dataVals2,lengthscale,initializationtype,parameterranges,use_fixedrnavariance);

fprintf(1,'createNdSimDisim step 1\n');

% normalization of POL2
dataVals1=dataVals1-min(dataVals1);
if mean(dataVals1.^2)>0,
  dataVals1=dataVals1/sqrt(mean(dataVals1.^2));
end;


% RNA variance estimate based on unnormalized read counts
rnavar_sigma=0.004541162;
rnavars = dataVals2+rnavar_sigma*dataVals2.^2;


% normalization of RNA and RNA variance
if ~isempty(dataVals2),
  if mean(dataVals2.^2)>0,
    normalizationfactor=sqrt(mean(dataVals2.^2));
    dataVals2=dataVals2/normalizationfactor;
    rnavars=rnavars/(normalizationfactor^2);
    fprintf('Normalizing!\n')
  end;
end;


fprintf(1,'createNdSimDisim step 2\n');


if ~isempty(dataVals2),
  dataValsMatrix=[dataVals1 dataVals2];
  dataVarsMatrix=[0*dataVals1 rnavars];
  numgenes=1;
else
  dataValsMatrix=[dataVals1];
  dataValsMatrix=[0*dataVals1];
  numgenes=0;
end;

if isempty(use_fixedrnavariance),
  use_fixedrnavariance=1;
end;



annotation=[];

options=struct();
options.includeNoise=1;
options.addPriors=0;
options.asynchronyType=1;
%options.optimiser='conjgrad';
%options.optimiser='scg';
options.optimiser='quasinew';
options.use_disimstartmean=1;
options.use_fixedrnavariance=use_fixedrnavariance;
options.fix=[];


if numgenes>0,
  inversewidth_index=1;
  pol2effectvar_index=2;
  rnadecay_index=3;
  rnaeffectvar_index=4;
  rnadelay_index=5;
  pol2noisevar_index=6;
  rnanoisevar_index=7;
  rnabasal_index=8;
  rnastartmean_index=9;
  pol2mean_index=10;
else
  inversewidth_index=1;
  pol2effectvar_index=2;
  pol2noisevar_index=3;
  rnadecay_index=nan;
  rnaeffectvar_index=nan;
  rnanoisevar_index=nan;
  rnabasal_index=nan;
  rnastartmean_index=nan;
  pol2mean_index=4;
end;

fprintf(1,'createNdSimDisim step 3\n');

if ~isempty(parameterranges),
  % use given parameter ranges
  inversewidth_range=parameterranges(1,:);
  pol2effectvar_range=parameterranges(2,:);
  pol2noisevar_range=parameterranges(3,:);
  rnadecay_range=parameterranges(4,:);
  rnaeffectvar_range=parameterranges(5,:);
  rnadelay_range=parameterranges(6,:);
  rnanoisevar_range=parameterranges(7,:);
  rnabasal_range=parameterranges(8,:);
  rnastartmean_range=parameterranges(9,:);
  pol2mean_range=parameterranges(10,:);  
else  
  % set parameter ranges
  inversewidth_range=1./([1280 5].^2);

  pol2effectvar_range=[0.0002 20]*var(dataVals1);
  if pol2effectvar_range(2)==pol2effectvar_range(1),
    pol2effectvar_range(2)=pol2effectvar_range(1)+1;
  end;
  
  pol2noisevar_range=[0.05 1]*var(dataVals1);
  if pol2noisevar_range(2)==pol2noisevar_range(1),
    pol2noisevar_range(2)=pol2noisevar_range(1)+1;
  end;
  
  rnadecay_range=[0.000001 1.5];
  rnaeffectvar_range=[0 20]; % *var(dataVals2)/var(dataVals1);
  %rnadelay_range=[0 800];
  rnadelay_range=[0 299];
  %rnadelay_range=[0 5];

  rnanoisevar_range=[0.05 1]*var(dataVals2);
  if rnanoisevar_range(2)==rnanoisevar_range(1),
    rnanoisevar_range(2)=rnanoisevar_range(1)+1;
  end;

  rnabasal_range=[0 10000];
  rnastartmean_range=[0 10000];
  pol2mean_range=[0 10000];
end;
  
transformsettings={};
if(~isnan(inversewidth_index)), transformsettings{inversewidth_index}=inversewidth_range; end;
if(~isnan(pol2effectvar_index)), transformsettings{pol2effectvar_index}=pol2effectvar_range; end;
if(~isnan(pol2noisevar_index)), transformsettings{pol2noisevar_index}=pol2noisevar_range; end;
if(~isnan(rnadecay_index)), transformsettings{rnadecay_index}=rnadecay_range; end;
if(~isnan(rnaeffectvar_index)), transformsettings{rnaeffectvar_index}=rnaeffectvar_range; end;
if(~isnan(rnadelay_index)), transformsettings{rnadelay_index}=rnadelay_range; end;
if(~isnan(rnanoisevar_index)), transformsettings{rnanoisevar_index}=rnanoisevar_range; end;
if(~isnan(rnabasal_index)), transformsettings{rnabasal_index}=rnabasal_range; end;
if(~isnan(rnastartmean_index)), transformsettings{rnastartmean_index}=rnastartmean_range; end;
if(~isnan(pol2mean_index)), transformsettings{pol2mean_index}=pol2mean_range; end;

transforminfo.indices=[inversewidth_index pol2effectvar_index rnadecay_index rnaeffectvar_index rnadelay_index pol2noisevar_index rnanoisevar_index rnabasal_index rnastartmean_index pol2mean_index];

transforminfo.settings={inversewidth_range,pol2effectvar_range,rnadecay_range,rnaeffectvar_range,rnadelay_range,pol2noisevar_range,rnanoisevar_range,rnabasal_range,rnastartmean_range,pol2mean_range};


gpsim3model = gpnddisimCreate(numgenes, 1, timevector, dataValsMatrix,dataVarsMatrix, options, annotation, transformsettings );


fprintf(1,'createNdSimDisim step 6\n');



if 1,

[pars,nams]=gpnddisimExtractParam(gpsim3model);



if(initializationtype==1),
  % initialization of inverse squared width
  if isempty(lengthscale),
    inversewidth_initval=1/(2*(15^2));
  else
    inversewidth_initval=1/(2*(lengthscale^2));
  end;
  pars(inversewidth_index)=sigmoidabTransform(inversewidth_initval, 'xtoa', inversewidth_range);
  
  % initialization of POL2 effect variance
  pol2effectvar_initval=0.04*var(dataVals1);
  pars(pol2effectvar_index)=sigmoidabTransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
  
  % initialization of POL2 noise variance
  pol2noisevar_initval=0.4*var(dataVals1);
  pars(pol2noisevar_index)=sigmoidabTransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

  % initialization of POL2 mean
  %pol2mean_initval=mean(dataVals1);
  I=find(timevector==min(timevector));    
  pol2mean_initval=mean(dataVals1(I));
  pars(pol2mean_index)=sigmoidabTransform(pol2mean_initval, 'xtoa', pol2mean_range);      
  
  if numgenes>0,
    % initialization of RNA decay
    rnadecay_initval=exp(-4);
    pars(rnadecay_index)=sigmoidabTransform(rnadecay_initval, 'xtoa', rnadecay_range);
    
    % initialization of RNA effect variance
    % idea: effect-kernel diagonal = RNAvar*POL2var*... = fraction of observed RNA variance
    %rnaeffectvar_initval=0.0002*var(dataVals2)/pol2noisevar_initval/(max(timevector)^3);
    rnaeffectvar_initval=0.03;
    S_initval=sqrt(rnaeffectvar_initval);
    pars(rnaeffectvar_index)=sigmoidabTransform(rnaeffectvar_initval, 'xtoa', rnaeffectvar_range);
    
    % initialization of RNA noise variance
    rnanoisevar_initval=0.4*var(dataVals2);
    pars(rnanoisevar_index)=sigmoidabTransform(rnanoisevar_initval, 'xtoa', rnanoisevar_range);
    
    % initialization of RNA basal rate
    % Idea: set (B+polmean*S)/D = RNAendvalue
    rnaendmean_initval=dataVals2(end);
    rnabasal_initval=rnaendmean_initval*rnadecay_initval-pol2mean_initval*S_initval;
    pars(rnabasal_index)=sigmoidabTransform(rnabasal_initval, 'xtoa', rnabasal_range);  
  
    % initialization of RNA start mean
    % model.disimStartMean(k)*exp(model.D(k)*(-min(predtimes))) = desiredval
    % --> disimStartMean = desiredval*exp(model.D(k)*(min(predtimes)))
    I=find(timevector==min(timevector));    
    rnastartmean_initval=mean(dataVals2(I));
    rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector));
    pars(rnastartmean_index)=sigmoidabTransform(rnastartmean_initval, 'xtoa', rnastartmean_range);
  
    % initialization of RNA delay
    rnadelay_initval=1e-3;
    pars(rnadelay_index)=sigmoidabTransform(rnadelay_initval, 'xtoa', rnadelay_range);
  end;
  
end;



if(initializationtype==2),
  % initialization of inverse squared width
  if isempty(lengthscale),
    inversewidth_initval=1/(2*(15^2));
  else
    inversewidth_initval=1/(2*(lengthscale^2));
  end;
  pars(inversewidth_index)=sigmoidabTransform(inversewidth_initval, 'xtoa', inversewidth_range);
  
  % initialization of POL2 effect variance
  pol2effectvar_initval=0.00015*var(dataVals1);
  pars(pol2effectvar_index)=sigmoidabTransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
  
  % initialization of POL2 noise variance
  pol2noisevar_initval=1.25*var(dataVals1);
  pars(pol2noisevar_index)=sigmoidabTransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

  % initialization of POL2 mean
  pol2mean_initval=mean(dataVals1);
  pars(pol2mean_index)=sigmoidabTransform(pol2mean_initval, 'xtoa', pol2mean_range);      
  
  if numgenes>0,
    % initialization of RNA decay
    rnadecay_initval=exp(-1);
    pars(rnadecay_index)=sigmoidabTransform(rnadecay_initval, 'xtoa', rnadecay_range);
    
    % initialization of RNA effect variance
    % idea: effect-kernel diagonal = RNAvar*POL2var*... = fraction of observed RNA variance    
    rnaeffectvar_initval=0.00000015*var(dataVals2)/pol2effectvar_initval
    S_initval=sqrt(rnaeffectvar_initval);
    pars(rnaeffectvar_index)=sigmoidabTransform(rnaeffectvar_initval, 'xtoa', rnaeffectvar_range)
        
    % initialization of RNA noise variance
    rnanoisevar_initval=1.25*var(dataVals2);
    pars(rnanoisevar_index)=sigmoidabTransform(rnanoisevar_initval, 'xtoa', rnanoisevar_range);
    
    % initialization of RNA basal rate
    % Idea: set (B+polmean*S)/D = RNAendvalue
    rnaendmean_initval=mean(dataVals2);
    rnabasal_initval=rnaendmean_initval*rnadecay_initval-pol2mean_initval*S_initval;
    pars(rnabasal_index)=sigmoidabTransform(rnabasal_initval, 'xtoa', rnabasal_range);  
%    rnadecay_initval
%    pol2mean_initval
%    S_initval

    % initialization of RNA start mean
    rnastartmean_initval=mean(dataVals2);
    rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector));
    pars(rnastartmean_index)=sigmoidabTransform(rnastartmean_initval, 'xtoa', rnastartmean_range);    

    % initialization of RNA delay
    rnadelay_initval=1e-3;
    pars(rnadelay_index)=sigmoidabTransform(rnadelay_initval, 'xtoa', rnadelay_range);        
  end;
end;




if(initializationtype==3),
  % initialization of inverse squared width
  if isempty(lengthscale),
    inversewidth_initval=1/(2*(15^2));
  else
    inversewidth_initval=1/(2*(lengthscale^2));
  end;
  pars(inversewidth_index)=sigmoidabTransform(inversewidth_initval, 'xtoa', inversewidth_range);
  
  % initialization of POL2 effect variance
  pol2effectvar_initval=0.1;
  pars(pol2effectvar_index)=sigmoidabTransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
  
  % initialization of POL2 noise variance
  pol2noisevar_initval=0.1;
  pars(pol2noisevar_index)=sigmoidabTransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

  % initialization of POL2 mean
  %pol2mean_initval=mean(dataVals1);
  pol2mean_initval=dataVals1(1);
  pars(pol2mean_index)=sigmoidabTransform(pol2mean_initval, 'xtoa', pol2mean_range);      
  
  if numgenes>0,
    % initialization of RNA decay
    rnadecay_initval=100;
    pars(rnadecay_index)=sigmoidabTransform(rnadecay_initval, 'xtoa', rnadecay_range);
    
    % initialization of RNA effect variance
    % idea: effect-kernel diagonal = RNAvar*POL2var*... = fraction of observed RNA variance
    rnaeffectvar_initval=100*100;
    S_initval=sqrt(rnaeffectvar_initval);
    pars(rnaeffectvar_index)=sigmoidabTransform(rnaeffectvar_initval, 'xtoa', rnaeffectvar_range);
    
    % initialization of RNA noise variance
    rnanoisevar_initval=0.1;
    pars(rnanoisevar_index)=sigmoidabTransform(rnanoisevar_initval, 'xtoa', rnanoisevar_range);
    
    % initialization of RNA basal rate
    % Idea: set (B+polmean*S)/D = RNAstartvalue
    rnabasal_initval=dataVals2(1)*rnadecay_initval - pol2mean_initval*S_initval;
    pars(rnabasal_index)=sigmoidabTransform(rnabasal_initval, 'xtoa', rnabasal_range);  
  
    % initialization of RNA start mean
    % idea: set to (B+polmean*S)/D
    rnastartmean_initval=(rnabasal_initval+pol2mean_initval*S_initval)/rnadecay_initval;
    rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector));
    pars(rnastartmean_index)=sigmoidabTransform(rnastartmean_initval, 'xtoa', rnastartmean_range);

    % initialization of RNA delay
    rnadelay_initval=1e-3;
    pars(rnadelay_index)=sigmoidabTransform(rnadelay_initval, 'xtoa', rnadelay_range);    
  end;
  
end;


if(initializationtype>3),
  bestpars=-1;
  bestll=-inf;

  ntrials=50;
  for itrial=1:ntrials,
    if isempty(lengthscale),
      % initialization of inverse squared width
      pars(inversewidth_index)=sigmoidabTransform(1/(2*(15^2)), 'xtoa', inversewidth_range);
    else
      pars(inversewidth_index)=sigmoidabTransform(1/(2*(lengthscale^2)), 'xtoa', inversewidth_range);  
    end;
    
    % initialization of POL2 effect variance
    pars(pol2effectvar_index)=sigmoidabTransform((rand^2)*var(dataVals1), 'xtoa', pol2effectvar_range);
    
    % initialization of POL2 noise variance
    pars(pol2noisevar_index)=sigmoidabTransform((rand^2)*var(dataVals1), 'xtoa', pol2noisevar_range);

    % initialization of POL2 mean
    pars(pol2mean_index)=sigmoidabTransform(rand, 'xtoa', pol2mean_range);
    
    if numgenes>0,
      % initialization of RNA decay
      rnadecay_initval=rnadecay_range(1)+(rand^4)*(rnadecay_range(2)-rnadecay_range(1));
      pars(rnadecay_index)=sigmoidabTransform(rnadecay_initval, 'xtoa', rnadecay_range);
      
      % initialization of RNA effect variance
      pars(rnaeffectvar_index)=sigmoidabTransform(rnaeffectvar_range(1)+(rand^4)*(rnaeffectvar_range(2)-rnaeffectvar_range(1)), 'xtoa', rnaeffectvar_range);
      
      % initialization of RNA noise variance
      pars(rnanoisevar_index)=sigmoidabTransform(rnanoisevar_range(1)+(rand^4)*(rnanoisevar_range(2)-rnanoisevar_range(1)), 'xtoa', rnanoisevar_range);
      
      % initialization of RNA basal rate
      pars(rnabasal_index)=sigmoidabTransform(rnabasal_range(1)+(rand^3)*(rnabasal_range(2)-rnabasal_range(1)), 'xtoa', rnabasal_range);  
      
      % initialization of RNA start mean
      rnastartmean_initval=rand;
      rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector));    
      pars(rnastartmean_index)=sigmoidabTransform(rnastartmean_initval, 'xtoa', rnastartmean_range);  

      % initialization of RNA delay
      % rnadelay_initval=rand*(max(timevector)-min(timevector))/2;
      rnadelay_initval=rnadelay_range(1)+(rand)*(rnadelay_range(2)-rnadelay_range(1))
      pars(rnadelay_index)=sigmoidabTransform(rnadelay_initval, 'xtoa', rnadelay_range);          
    end;    
	
    gpsim3model=gpnddisimExpandParam(gpsim3model,pars);
    templl=gpnddisimLogLikelihood(gpsim3model)
    if templl>bestll,
      bestpars=pars;
      bestll=templl;
    end;
  end;
  pars=bestpars;  
end;



%nams
%pars
%pause
gpsim3model=gpnddisimExpandParam(gpsim3model,pars);

%fprintf(1,'Trying these parameters:\n');
[pars,name]=gpnddisimExtractParam(gpsim3model);
%nams
%pars

if 0,
  transformedpars=pars;
  if(~isnan(pol2decay_index)), 
    tempindex=pol2decay_index;
    temprange=pol2decay_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(inversewidth_index)), 
    tempindex=inversewidth_index;
    temprange=inversewidth_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(pol2effectvar_index)), 
    tempindex=pol2effectvar_index;
    temprange=pol2effectvar_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(pol2noisevar_index)), 
    tempindex=pol2noisevar_index;
    temprange=pol2noisevar_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnadecay_index)), 
    tempindex=rnadecay_index;
    temprange=rnadecay_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnaeffectvar_index)), 
    tempindex=rnaeffectvar_index;
    temprange=rnaeffectvar_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rbfvar_index)), 
    tempindex=rbfvar_index;
    temprange=rbfvar_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnanoisevar_index)), 
    tempindex=rnanoisevar_index;
    temprange=rnanoisevar_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnabasal_index)), 
    tempindex=rnabasal_index;
    temprange=rnabasal_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnastartmean_index)), 
    tempindex=rnastartmean_index;
    temprange=rnastartmean_range;
    transformedpars(tempindex)=sigmoidabTransform(pars(tempindex),'atox',temprange); 
  end;

  %transformedpars
end;

end;

fprintf(1,'createSimDisim done\n');
