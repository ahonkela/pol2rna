function [gpsim3model,transforminfo]= ...
    createNdSimDisim_celltimes_newdata3(timevector,dataVals,lengthscale,initializationtype,parameterranges,use_fixedrnavariance,addPriors,rnavars,donotnormalize,uniformPriors);

DEBUG=0;

if DEBUG,
  fprintf(1,'createNdSimDisim step 1\n');
end

% normalization of POL2
dataVals1=dataVals{1};
if ~donotnormalize,
  if isempty(dataVals1)==0,
    dataVals1=dataVals1-min(dataVals1);
    if mean(dataVals1.^2)>0,
      dataVals1=dataVals1/sqrt(mean(dataVals1.^2));
    end;
  end;
end


if length(dataVals)>1,
  dataVals2=dataVals{2};
else
  dataVals2=[];
end;

% if isempty(dataVals2)==0,
%   % RNA variance estimate based on unnormalized read counts
%   rnavar_sigma=0.004238473;
%   tempvals = max(1,dataVals2);
%   rnavars = tempvals+rnavar_sigma*tempvals.^2;
% else
%   rnavar_sigma=0.004238473;
%   rnavars = [];
% end;

% normalization of RNA and RNA variance
if ~donotnormalize && ~isempty(dataVals2),
  if mean(dataVals2.^2)>0,
    normalizationfactor=sqrt(mean(dataVals2.^2));
    dataVals2=dataVals2/normalizationfactor;
    if iscell(rnavars),
      rnavars{2}=rnavars{2}/(normalizationfactor^2);
    else
      rnavars=rnavars/(normalizationfactor^2);
    end
    fprintf('Normalizing!\n')
  end;
end;

if DEBUG,
  fprintf(1,'createNdSimDisim step 2\n');
end


%if ~isempty(dataVals2),
if length(dataVals)>1,
  % dataValsMatrix=[dataVals1 dataVals2];
  % dataVarsMatrix=[0*dataVals1 rnavars];
  dataValsCell={dataVals1, dataVals2};
  if iscell(rnavars),
    dataVarsCell=rnavars;
  else
    dataVarsCell={0*dataVals1, rnavars};
  end
  numgenes=1;
else
  % dataValsMatrix=[dataVals1];
  % dataVarsMatrix=[0*dataVals1];
  dataValsCell={dataVals1};
  if iscell(rnavars),
    dataVarsCell=rnavars;
  else
    dataVarsCell={0*dataVals1};
  end
  numgenes=0;
end;

if isempty(use_fixedrnavariance),
  use_fixedrnavariance=1;
end;



annotation=[];

options=struct();
options.includeNoise=1;
% Remember to add priors!
if ~isempty(addPriors) && addPriors,
  options.addPriors=1;
else
  options.addPriors=0;
end
options.asynchronyType=1;
%options.optimiser='conjgrad';
%options.optimiser='scg';
options.optimiser='quasinew';
options.use_disimstartmean=1;
options.use_fixedrnavariance=use_fixedrnavariance;
options.fix=[];
options.uniformPriors = uniformPriors;

if uniformPriors,
  mytransform = str2func('identityTransform');
else
  mytransform = str2func('sigmoidabTransform');
end

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

if DEBUG,
  fprintf(1,'createNdSimDisim step 3\n');
end

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

  if isempty(dataVals1)==0,
    pol2effectvar_range=[0.0002 1]*var(dataVals1);
  else
    pol2effectvar_range=[0.0002 1];
  end;
  if pol2effectvar_range(2)==pol2effectvar_range(1),
    pol2effectvar_range(2)=pol2effectvar_range(1)+1;
  end;
  
  if isempty(dataVals1)==0,
    pol2noisevar_range=[0.05 1]*var(dataVals1);
  else
    pol2noisevar_range=[0.05 1];
  end;
  if pol2noisevar_range(2)==pol2noisevar_range(1),
    pol2noisevar_range(2)=pol2noisevar_range(1)+1;
  end;
  
  rnadecay_range=[0.000001 log(2)];
  rnaeffectvar_range=[1e-6 1]; % *var(dataVals2)/var(dataVals1);
  %rnadelay_range=[0 800];
  rnadelay_range=[0 299];
  %rnadelay_range=[0 5];

  if isempty(dataVals2)==0,
    rnanoisevar_range=[0.05 1]*var(dataVals2);
  else
    rnanoisevar_range=[0.05 1];
  end;
  if rnanoisevar_range(2)==rnanoisevar_range(1),
    rnanoisevar_range(2)=rnanoisevar_range(1)+1;
  end;

  rnabasal_range=[0 1];
  rnastartmean_range=[0 2];
  pol2mean_range=[0 1];
end;

if DEBUG,
  fprintf(1,'createNdSimDisim step 4\n');
end
  
transformsettings={};
if(~isnan(inversewidth_index)), transformsettings{inversewidth_index}=inversewidth_range; end;
if(~isnan(pol2effectvar_index)), transformsettings{pol2effectvar_index}=pol2effectvar_range; end;
if(~isnan(pol2noisevar_index)), transformsettings{pol2noisevar_index}=pol2noisevar_range; end;
if(~isnan(rnadecay_index)), transformsettings{rnadecay_index}=rnadecay_range; end;
if(~isnan(rnaeffectvar_index)), transformsettings{rnaeffectvar_index}=rnaeffectvar_range; end;

if DEBUG,
  fprintf(1,'createNdSimDisim step 4a\n');
  numgenes
  rnadelay_index
  rnadelay_range
end

if(~isnan(rnadelay_index)), transformsettings{rnadelay_index}=rnadelay_range; end;

if DEBUG,
  fprintf(1,'createNdSimDisim step 4a1\n');
end

if(~isnan(rnanoisevar_index)), transformsettings{rnanoisevar_index}=rnanoisevar_range; end;
if(~isnan(rnabasal_index)), transformsettings{rnabasal_index}=rnabasal_range; end;

if DEBUG,
  fprintf(1,'createNdSimDisim step 4a2\n');
end

if(~isnan(rnastartmean_index)), transformsettings{rnastartmean_index}=rnastartmean_range; end;
if(~isnan(pol2mean_index)), transformsettings{pol2mean_index}=pol2mean_range; end;

if DEBUG,
  fprintf(1,'createNdSimDisim step 4b\n');
end

transforminfo.indices=[inversewidth_index pol2effectvar_index rnadecay_index rnaeffectvar_index rnadelay_index pol2noisevar_index rnanoisevar_index rnabasal_index rnastartmean_index pol2mean_index];

transforminfo.settings={inversewidth_range,pol2effectvar_range,rnadecay_range,rnaeffectvar_range,rnadelay_range,pol2noisevar_range,rnanoisevar_range,rnabasal_range,rnastartmean_range,pol2mean_range};

if DEBUG,
  fprintf(1,'createNdSimDisim step 5\n');
end

gpsim3model = gpnddisimCreate(numgenes, 1, timevector, dataValsCell,dataVarsCell, options, annotation, transformsettings );


if DEBUG
  fprintf(1,'createNdSimDisim step 6\n');
end


if 1,

[pars,nams]=gpnddisimExtractParam(gpsim3model);

%pars
%nams
%pause



if(initializationtype==1),
  % initialization of inverse squared width
  if isempty(lengthscale),
    inversewidth_initval=1/(2*(15^2));
  else
    inversewidth_initval=1/(2*(lengthscale^2));
  end;
  pars(inversewidth_index)=mytransform(inversewidth_initval, 'xtoa', inversewidth_range);
  
  % initialization of POL2 effect variance
  if isempty(dataVals1)==0,
    pol2effectvar_initval=0.04*var(dataVals1);
  else
    pol2effectvar_initval=0.04;
  end;
  pars(pol2effectvar_index)=mytransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
  
  % initialization of POL2 noise variance
  if isempty(dataVals1)==0,
    pol2noisevar_initval=0.4*var(dataVals1);
  else
    pol2noisevar_initval=0.4;
  end;
  pars(pol2noisevar_index)=mytransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

  % initialization of POL2 mean
  %pol2mean_initval=mean(dataVals1);
  if isempty(dataVals1)==0,
    I=find(timevector{1}==min(timevector{1}));
    pol2mean_initval=mean(dataVals1(I));
  else
   pol2mean_initval=0.1;
  end;
  pars(pol2mean_index)=mytransform(pol2mean_initval, 'xtoa', pol2mean_range);      
  
  if numgenes>0,
    % initialization of RNA decay
    rnadecay_initval=exp(-4);
    pars(rnadecay_index)=mytransform(rnadecay_initval, 'xtoa', rnadecay_range);
    
    % initialization of RNA effect variance
    % idea: effect-kernel diagonal = RNAvar*POL2var*... = fraction of observed RNA variance
    %rnaeffectvar_initval=0.0002*var(dataVals2)/pol2noisevar_initval/(max(timevector)^3);
    rnaeffectvar_initval=0.03;
    S_initval=sqrt(rnaeffectvar_initval);
    pars(rnaeffectvar_index)=mytransform(rnaeffectvar_initval, 'xtoa', rnaeffectvar_range);
    
    % initialization of RNA noise variance
    if isempty(dataVals2)==0,
      rnanoisevar_initval=0.4*var(dataVals2);
    else
      rnanoisevar_initval=0.4;
    end;
    pars(rnanoisevar_index)=mytransform(rnanoisevar_initval, 'xtoa', rnanoisevar_range);
    
    % initialization of RNA basal rate
    % Idea: set (B+polmean*S)/D = RNAendvalue
    if isempty(dataVals2)==0, 
      rnaendmean_initval=dataVals2(end);
    else
      rnaendmean_initval=0.1;
    end;
    rnabasal_initval=rnaendmean_initval*rnadecay_initval-pol2mean_initval*S_initval;
    pars(rnabasal_index)=mytransform(rnabasal_initval, 'xtoa', rnabasal_range);  
  
    % initialization of RNA start mean
    % model.disimStartMean(k)*exp(model.D(k)*(-min(predtimes))) = desiredval
    % --> disimStartMean = desiredval*exp(model.D(k)*(min(predtimes)))
    if isempty(dataVals2)==0,
      I=find(timevector{2}==min(timevector{2}));    
      rnastartmean_initval=mean(dataVals2(I));
      rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector{2}));
    else
      rnastartmean_initval=0.1;
    end;
    pars(rnastartmean_index)=mytransform(rnastartmean_initval, 'xtoa', rnastartmean_range);
  
    % initialization of RNA delay
    rnadelay_initval=1e-3;
    pars(rnadelay_index)=mytransform(rnadelay_initval, 'xtoa', rnadelay_range);
  end;
  
end;



if(initializationtype==2),
  % initialization of inverse squared width
  if isempty(lengthscale),
    inversewidth_initval=1/(2*(15^2));
  else
    inversewidth_initval=1/(2*(lengthscale^2));
  end;
  pars(inversewidth_index)=mytransform(inversewidth_initval, 'xtoa', inversewidth_range);
  
  % initialization of POL2 effect variance
  if isempty(dataVals1)==0,
    pol2effectvar_initval=0.00015*var(dataVals1);
  else
    pol2effectvar_initval=0.00015;
  end;
  pars(pol2effectvar_index)=mytransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
  
  % initialization of POL2 noise variance
  if isempty(dataVals1)==0,
    pol2noisevar_initval=1.25*var(dataVals1);
  else
    pol2noisevar_initval=1.25;
  end;
  pars(pol2noisevar_index)=mytransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

  % initialization of POL2 mean
  if isempty(dataVals1)==0,
    pol2mean_initval=mean(dataVals1);
  else
    pol2mean_initval=0.1;
  end;
  pars(pol2mean_index)=mytransform(pol2mean_initval, 'xtoa', pol2mean_range);      
  
  if numgenes>0,
    % initialization of RNA decay
    rnadecay_initval=exp(-1);
    pars(rnadecay_index)=mytransform(rnadecay_initval, 'xtoa', rnadecay_range);
    
    % initialization of RNA effect variance
    % idea: effect-kernel diagonal = RNAvar*POL2var*... = fraction of observed RNA variance    
    if isempty(dataVals2)==0,
      rnaeffectvar_initval=0.00000015*var(dataVals2)/pol2effectvar_initval;
    else
      rnaeffectvar_initval=0.00000015/pol2effectvar_initval;
    end;
    S_initval=sqrt(rnaeffectvar_initval);
    pars(rnaeffectvar_index)=mytransform(rnaeffectvar_initval, 'xtoa', rnaeffectvar_range)
        
    % initialization of RNA noise variance
    if isempty(dataVals2)==0,
      rnanoisevar_initval=1.25*var(dataVals2);
    else
      rnanoisevar_initval=1.25;
    end;
    pars(rnanoisevar_index)=mytransform(rnanoisevar_initval, 'xtoa', rnanoisevar_range);
    
    % initialization of RNA basal rate
    % Idea: set (B+polmean*S)/D = RNAendvalue
    if isempty(dataVals2)==0,
      rnaendmean_initval=mean(dataVals2);
    else
      rnaendmean_initval=0.1;
    end;
    rnabasal_initval=rnaendmean_initval*rnadecay_initval-pol2mean_initval*S_initval;
    pars(rnabasal_index)=mytransform(rnabasal_initval, 'xtoa', rnabasal_range);  
%    rnadecay_initval
%    pol2mean_initval
%    S_initval

    % initialization of RNA start mean
    if isempty(dataVals2)==0 && 0,
      rnastartmean_initval=mean(dataVals2);
      rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector{2}));
    else
      rnastartmean_initval=0.1;
    end;
    pars(rnastartmean_index)=mytransform(rnastartmean_initval, 'xtoa', rnastartmean_range);    

    % initialization of RNA delay
    rnadelay_initval=1e-3;
    pars(rnadelay_index)=mytransform(rnadelay_initval, 'xtoa', rnadelay_range);        
  end;
end;




if(initializationtype==3),
  % initialization of inverse squared width
  if isempty(lengthscale),
    inversewidth_initval=1/(2*(15^2));
  else
    inversewidth_initval=1/(2*(lengthscale^2));
  end;
  pars(inversewidth_index)=mytransform(inversewidth_initval, 'xtoa', inversewidth_range);
  
  % initialization of POL2 effect variance
  pol2effectvar_initval=0.1;
  pars(pol2effectvar_index)=mytransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
  
  % initialization of POL2 noise variance
  pol2noisevar_initval=0.1;
  pars(pol2noisevar_index)=mytransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

  % initialization of POL2 mean
  %pol2mean_initval=mean(dataVals1);
  if isempty(dataVals1)==0,
    pol2mean_initval=dataVals1(1);
  else
    pol2mean_initval=0.1;
  end;
  pars(pol2mean_index)=mytransform(pol2mean_initval, 'xtoa', pol2mean_range);      
  
  if numgenes>0,
    % initialization of RNA decay
    rnadecay_initval=exp(-2);
    pars(rnadecay_index)=mytransform(rnadecay_initval, 'xtoa', rnadecay_range);
    
    % initialization of RNA effect variance
    % idea: effect-kernel diagonal = RNAvar*POL2var*... = fraction of observed RNA variance
    rnaeffectvar_initval=100*100;
    S_initval=sqrt(rnaeffectvar_initval);
    pars(rnaeffectvar_index)=mytransform(rnaeffectvar_initval, 'xtoa', rnaeffectvar_range);
    
    % initialization of RNA noise variance
    rnanoisevar_initval=0.1;
    pars(rnanoisevar_index)=mytransform(rnanoisevar_initval, 'xtoa', rnanoisevar_range);
    
    % initialization of RNA basal rate
    % Idea: set (B+polmean*S)/D = RNAstartvalue
    if isempty(dataVals2)==0,
      rnabasal_initval=dataVals2(1)*rnadecay_initval - pol2mean_initval*S_initval;
    else
      rnabasal_initval=0.1*rnadecay_initval - pol2mean_initval*S_initval;
    end;
    pars(rnabasal_index)=mytransform(rnabasal_initval, 'xtoa', rnabasal_range);  
  
    % initialization of RNA start mean
    % idea: set to (B+polmean*S)/D
    if isempty(dataVals2)==0,
      rnastartmean_initval=(rnabasal_initval+pol2mean_initval*S_initval)/rnadecay_initval;
      rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector{2}));
    else
      rnastartmean_initval=(rnabasal_initval+pol2mean_initval*S_initval)/rnadecay_initval;
      rnastartmean_initval=rnastartmean_initval;
    end;
    pars(rnastartmean_index)=mytransform(rnastartmean_initval, 'xtoa', rnastartmean_range);

    % initialization of RNA delay
    rnadelay_initval=1e-3;
    pars(rnadelay_index)=mytransform(rnadelay_initval, 'xtoa', rnadelay_range);    
  end;
  
end;




if(initializationtype==4),
  % initialization of inverse squared width
  if isempty(lengthscale),
    inversewidth_initval=1/(2*(15^2));
  else
    inversewidth_initval=1/(2*(lengthscale^2));
  end;
  pars(inversewidth_index)=mytransform(inversewidth_initval, 'xtoa', inversewidth_range);
  
  % initialization of POL2 effect variance
  if isempty(dataVals1)==0,
    pol2effectvar_initval=0.01*var(dataVals1);
  else
    pol2effectvar_initval=0.1;
  end;
  pars(pol2effectvar_index)=mytransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
  
  % initialization of POL2 noise variance
  if isempty(dataVals1)==0,
    pol2noisevar_initval=0.01*var(dataVals1);
  else
    pol2noisevar_initval=0.01;
  end;
  pars(pol2noisevar_index)=mytransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

  % initialization of POL2 mean
  %pol2mean_initval=mean(dataVals1);
  if isempty(dataVals1)==0,
    I=find(timevector{1}==min(timevector{1}));
    pol2mean_initval=mean(dataVals1(I));
  else
   pol2mean_initval=0.1;
  end;
  pars(pol2mean_index)=mytransform(pol2mean_initval, 'xtoa', pol2mean_range);      
  
  if numgenes>0,
    % initialization of RNA decay
    rnadecay_initval=exp(-4);
    pars(rnadecay_index)=mytransform(rnadecay_initval, 'xtoa', rnadecay_range);
    
    % initialization of RNA effect variance
    % idea: effect-kernel diagonal = RNAvar*POL2var*... = fraction of observed RNA variance
    %rnaeffectvar_initval=0.0002*var(dataVals2)/pol2noisevar_initval/(max(timevector)^3);
    rnaeffectvar_initval=0.0001;
    S_initval=sqrt(rnaeffectvar_initval);
    pars(rnaeffectvar_index)=mytransform(rnaeffectvar_initval, 'xtoa', rnaeffectvar_range);
    
    % initialization of RNA noise variance
    if isempty(dataVals2)==0,
      rnanoisevar_initval=0.01*var(dataVals2);
    else
      rnanoisevar_initval=0.01;
    end;
    pars(rnanoisevar_index)=mytransform(rnanoisevar_initval, 'xtoa', rnanoisevar_range);
    
    % initialization of RNA basal rate
    % Idea: set (B+polmean*S)/D = RNAendvalue
    if isempty(dataVals2)==0, 
      rnaendmean_initval=dataVals2(end);
    else
      rnaendmean_initval=0.1;
    end;
    rnabasal_initval=rnaendmean_initval*rnadecay_initval-pol2mean_initval*S_initval;
    pars(rnabasal_index)=mytransform(rnabasal_initval, 'xtoa', rnabasal_range);  
  
    % initialization of RNA start mean
    % model.disimStartMean(k)*exp(model.D(k)*(-min(predtimes))) = desiredval
    % --> disimStartMean = desiredval*exp(model.D(k)*(min(predtimes)))
    if isempty(dataVals2)==0,
      I=find(timevector{2}==min(timevector{2}));    
      rnastartmean_initval=mean(dataVals2(I));
      rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector{2}));
    else
      rnastartmean_initval=0.1;
    end;
    pars(rnastartmean_index)=mytransform(rnastartmean_initval, 'xtoa', rnastartmean_range);
  
    % initialization of RNA delay
    rnadelay_initval=1e-3;
    pars(rnadelay_index)=mytransform(rnadelay_initval, 'xtoa', rnadelay_range);
  end;
  
end;




if(initializationtype==5),
  % initialization of inverse squared width
  if isempty(lengthscale),
    inversewidth_initval=1/(2*(15^2));
  else
    inversewidth_initval=1/(2*(lengthscale^2));
  end;
  pars(inversewidth_index)=mytransform(inversewidth_initval, 'xtoa', inversewidth_range);
  
  % initialization of POL2 effect variance
  if isempty(dataVals1)==0,
    pol2effectvar_initval=0.001*var(dataVals1);
  else
    pol2effectvar_initval=0.01;
  end;
  pars(pol2effectvar_index)=mytransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
  
  % initialization of POL2 noise variance
  if isempty(dataVals1)==0,
    pol2noisevar_initval=0.9*var(dataVals1);
  else
    pol2noisevar_initval=0.9;
  end;
  pars(pol2noisevar_index)=mytransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

  % initialization of POL2 mean
  %pol2mean_initval=mean(dataVals1);
  if isempty(dataVals1)==0,
    I=find(timevector{1}==min(timevector{1}));
    pol2mean_initval=mean(dataVals1(I));
  else
   pol2mean_initval=0.1;
  end;
  pars(pol2mean_index)=mytransform(pol2mean_initval, 'xtoa', pol2mean_range);      
  
  if numgenes>0,
    % initialization of RNA decay
    rnadecay_initval=exp(-4);
    pars(rnadecay_index)=mytransform(rnadecay_initval, 'xtoa', rnadecay_range);
    
    % initialization of RNA effect variance
    % idea: effect-kernel diagonal = RNAvar*POL2var*... = fraction of observed RNA variance
    %rnaeffectvar_initval=0.0002*var(dataVals2)/pol2noisevar_initval/(max(timevector)^3);
    rnaeffectvar_initval=0.00001;
    S_initval=sqrt(rnaeffectvar_initval);
    pars(rnaeffectvar_index)=mytransform(rnaeffectvar_initval, 'xtoa', rnaeffectvar_range);
    
    % initialization of RNA noise variance
    if isempty(dataVals2)==0,
      rnanoisevar_initval=0.9*var(dataVals2);
    else
      rnanoisevar_initval=0.9;
    end;
    pars(rnanoisevar_index)=mytransform(rnanoisevar_initval, 'xtoa', rnanoisevar_range);
    
    % initialization of RNA basal rate
    % Idea: set (B+polmean*S)/D = RNAendvalue
    if isempty(dataVals2)==0, 
      rnaendmean_initval=dataVals2(end);
    else
      rnaendmean_initval=0.1;
    end;
    rnabasal_initval=rnaendmean_initval*rnadecay_initval-pol2mean_initval*S_initval;
    pars(rnabasal_index)=mytransform(rnabasal_initval, 'xtoa', rnabasal_range);  
  
    % initialization of RNA start mean
    % model.disimStartMean(k)*exp(model.D(k)*(-min(predtimes))) = desiredval
    % --> disimStartMean = desiredval*exp(model.D(k)*(min(predtimes)))
    if isempty(dataVals2)==0,
      I=find(timevector{2}==min(timevector{2}));    
      rnastartmean_initval=mean(dataVals2(I));
      rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector{2}));
    else
      rnastartmean_initval=0.1;
    end;
    pars(rnastartmean_index)=mytransform(rnastartmean_initval, 'xtoa', rnastartmean_range);
  
    % initialization of RNA delay
    rnadelay_initval=1e-3;
    pars(rnadelay_index)=mytransform(rnadelay_initval, 'xtoa', rnadelay_range);
  end;
  
end;


if any(isnan(pars)),
    warning('NaN parameters');
    pars(isnan(pars)) = 0;
end


if(initializationtype==6),
  bestpars=-1;
  bestll=-inf;

  ntrials=50;
  for itrial=1:ntrials,
    if isempty(lengthscale),
      % initialization of inverse squared width
      pars(inversewidth_index)=mytransform(1/(2*(15^2)), 'xtoa', inversewidth_range);
    else
      pars(inversewidth_index)=mytransform(1/(2*(lengthscale^2)), 'xtoa', inversewidth_range);  
    end;
    
    % initialization of POL2 effect variance
    if isempty(dataVals1)==0,
      pol2effectvar_initval=(rand^2)*var(dataVals1);
    else
      pol2effectvar_initval=(rand^2);
    end;
    pars(pol2effectvar_index)=mytransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
    
    % initialization of POL2 noise variance
    if isempty(dataVals1)==0,
      pol2noisevar_initval=(rand^2)*var(dataVals1);
    else
      pol2noisevar_initval=(rand^2);
    end;
    pars(pol2noisevar_index)=mytransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

    % initialization of POL2 mean
    pars(pol2mean_index)=mytransform(rand, 'xtoa', pol2mean_range);
    
    if numgenes>0,
      % initialization of RNA decay
      rnadecay_initval=rnadecay_range(1)+(rand^4)*(rnadecay_range(2)-rnadecay_range(1));
      pars(rnadecay_index)=mytransform(rnadecay_initval, 'xtoa', rnadecay_range);
      
      % initialization of RNA effect variance
      pars(rnaeffectvar_index)=mytransform(rnaeffectvar_range(1)+(rand^4)*(rnaeffectvar_range(2)-rnaeffectvar_range(1)), 'xtoa', rnaeffectvar_range);
      
      % initialization of RNA noise variance
      pars(rnanoisevar_index)=mytransform(rnanoisevar_range(1)+(rand^4)*(rnanoisevar_range(2)-rnanoisevar_range(1)), 'xtoa', rnanoisevar_range);
      
      % initialization of RNA basal rate
      pars(rnabasal_index)=mytransform(rnabasal_range(1)+(rand^3)*(rnabasal_range(2)-rnabasal_range(1)), 'xtoa', rnabasal_range);  
      
      % initialization of RNA start mean
      if isempty(dataVals2)==0,
        rnastartmean_initval=rand;
        rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector{2}));    
      else
        rnastartmean_initval=rand;
      end;
      pars(rnastartmean_index)=mytransform(rnastartmean_initval, 'xtoa', rnastartmean_range);  

      % initialization of RNA delay
      % rnadelay_initval=rand*(max(timevector)-min(timevector))/2;
      rnadelay_initval=rnadelay_range(1)+(rand)*(rnadelay_range(2)-rnadelay_range(1));
      pars(rnadelay_index)=mytransform(rnadelay_initval, 'xtoa', rnadelay_range);          
    end;    
	
    gpsim3model=gpnddisimExpandParam(gpsim3model,pars);
    templl=gpnddisimLogLikelihood(gpsim3model);
    if templl>bestll,
      bestpars=pars;
      bestll=templl;
    end;
  end;
  pars=bestpars;  
end;


if(initializationtype>6),
  rng(initializationtype);
  bestpars=-1;
  bestll=-inf;

  ntrials=3;
  for itrial=1:ntrials,
    if uniformPriors,
      pars0 = gpnddisimExtractParam(gpsim3model);
      pars = rand(size(pars0));
      for k=1:length(pars),
        pars(k) = transformsettings{k}(1) + ...
                  pars(k) * (transformsettings{k}(2)-transformsettings{k}(1));
      end
    else
    if isempty(lengthscale),
      % initialization of inverse squared width
      pars(inversewidth_index)=mytransform(1/(2*(15^2)), 'xtoa', inversewidth_range);
    else
      pars(inversewidth_index)=mytransform(1/(2*(lengthscale^2)), 'xtoa', inversewidth_range);  
    end;
    
    % initialization of POL2 effect variance
    if isempty(dataVals1)==0,
      pol2effectvar_initval=(rand^2)*var(dataVals1);
    else
      pol2effectvar_initval=(rand^2);
    end;
    pars(pol2effectvar_index)=mytransform(pol2effectvar_initval, 'xtoa', pol2effectvar_range);
    
    % initialization of POL2 noise variance
    if isempty(dataVals1)==0,
      pol2noisevar_initval=(rand^2)*var(dataVals1);
    else
      pol2noisevar_initval=(rand^2);
    end;
    pars(pol2noisevar_index)=mytransform(pol2noisevar_initval, 'xtoa', pol2noisevar_range);

    % initialization of POL2 mean
    pars(pol2mean_index)=mytransform(rand, 'xtoa', pol2mean_range);
    
    if numgenes>0,
      % initialization of RNA decay
      rnadecay_initval=rnadecay_range(1)+(rand^4)*(rnadecay_range(2)-rnadecay_range(1));
      pars(rnadecay_index)=mytransform(rnadecay_initval, 'xtoa', rnadecay_range);
      
      % initialization of RNA effect variance
      pars(rnaeffectvar_index)=mytransform(rnaeffectvar_range(1)+(rand^4)*(rnaeffectvar_range(2)-rnaeffectvar_range(1)), 'xtoa', rnaeffectvar_range);
      
      % initialization of RNA noise variance
      pars(rnanoisevar_index)=mytransform(rnanoisevar_range(1)+(rand^4)*(rnanoisevar_range(2)-rnanoisevar_range(1)), 'xtoa', rnanoisevar_range);
      
      % initialization of RNA basal rate
      pars(rnabasal_index)=mytransform(rnabasal_range(1)+(rand^3)*(rnabasal_range(2)-rnabasal_range(1)), 'xtoa', rnabasal_range);  
      
      % initialization of RNA start mean
      if isempty(dataVals2)==0,
        rnastartmean_initval=rand;
        rnastartmean_initval=rnastartmean_initval*exp(rnadecay_initval*min(timevector{2}));    
      else
        rnastartmean_initval=rand;
      end;
      pars(rnastartmean_index)=mytransform(rnastartmean_initval, 'xtoa', rnastartmean_range);  

      % initialization of RNA delay
      % rnadelay_initval=rand*(max(timevector)-min(timevector))/2;
      rnadelay_initval=rnadelay_range(1)+(rand^3)*(rnadelay_range(2)-rnadelay_range(1));
      pars(rnadelay_index)=mytransform(rnadelay_initval, 'xtoa', rnadelay_range);          
    end;    
    end   % if ~uniformPriors
	
    gpsim3model=gpnddisimExpandParam(gpsim3model,pars);
    templl=gpnddisimLogLikelihood(gpsim3model);
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
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(inversewidth_index)), 
    tempindex=inversewidth_index;
    temprange=inversewidth_range;
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(pol2effectvar_index)), 
    tempindex=pol2effectvar_index;
    temprange=pol2effectvar_range;
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(pol2noisevar_index)), 
    tempindex=pol2noisevar_index;
    temprange=pol2noisevar_range;
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnadecay_index)), 
    tempindex=rnadecay_index;
    temprange=rnadecay_range;
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnaeffectvar_index)), 
    tempindex=rnaeffectvar_index;
    temprange=rnaeffectvar_range;
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rbfvar_index)), 
    tempindex=rbfvar_index;
    temprange=rbfvar_range;
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnanoisevar_index)), 
    tempindex=rnanoisevar_index;
    temprange=rnanoisevar_range;
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnabasal_index)), 
    tempindex=rnabasal_index;
    temprange=rnabasal_range;
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;
  if(~isnan(rnastartmean_index)), 
    tempindex=rnastartmean_index;
    temprange=rnastartmean_range;
    transformedpars(tempindex)=mytransform(pars(tempindex),'atox',temprange); 
  end;

  %transformedpars
end;

end;

if DEBUG,
  fprintf(1,'createSimDisim done\n');
end
