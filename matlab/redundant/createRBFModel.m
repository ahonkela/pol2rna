function [pol2model,transforminfo] = createRBFModel(timevector,dataVals1,lengthscale,initializationtype);

% normalization of POL2
%dataVals1=dataVals1-dataVals1(1);
dataVals1=dataVals1-mean(dataVals1);
dataVals1=dataVals1/sqrt(var(dataVals1));


% set parameter ranges
inversewidth_index=1;
pol2effectvar_index=2;
pol2noisevar_index=3;
inversewidth_range=1./(2*([320 5].^2));
pol2effectvar_range=[0 10]*var(dataVals1);
pol2noisevar_range=[0.02 2]*var(dataVals1);

transformsettings={};
if(~isnan(inversewidth_index)), transformsettings{inversewidth_index}=inversewidth_range; end;
if(~isnan(pol2effectvar_index)), transformsettings{pol2effectvar_index}=pol2effectvar_range; end;
if(~isnan(pol2noisevar_index)), transformsettings{pol2noisevar_index}=pol2noisevar_range; end;

transforminfo.indices=[inversewidth_index pol2effectvar_index pol2noisevar_index];
transforminfo.settings={inversewidth_range,pol2effectvar_range,pol2noisevar_range};





options=struct();
options.includeNoise=1;
options.addPriors=1;
%options.optimiser='conjgrad';
options.optimiser='scg';
options.fix=[];
annotation=[];
pol2model = gpasimTempRBFCreate(0, 1, timevector, dataVals1, 0*dataVals1, options, annotation, transformsettings );
[pars,nams]=gpdisimExtractParam(pol2model);

%inversewidth_index=1;
%pol2effectvar_index=2;
%pol2noisevar_index=3;

if (initializationtype==1),
  if isempty(lengthscale),
    % initialization of inverse squared width
    pars(inversewidth_index)=sigmoidabTransform(1/(2*(15^2)), 'xtoa', inversewidth_range);
  else
    pars(inversewidth_index)=sigmoidabTransform(1/(2*(lengthscale^2)), 'xtoa', inversewidth_range);  
  end;
  
  % initialization of POL2 effect variance
  pars(pol2effectvar_index)=sigmoidabTransform(0.2*var(dataVals1), 'xtoa', pol2effectvar_range);
  % initialization of POL2 noise variance
  pars(pol2noisevar_index)=sigmoidabTransform(0.5*var(dataVals1), 'xtoa', pol2noisevar_range);
end;

if (initializationtype==2),
  if isempty(lengthscale),
    % initialization of inverse squared width
    pars(inversewidth_index)=sigmoidabTransform(1/(2*(15^2)), 'xtoa', inversewidth_range);
  else
    pars(inversewidth_index)=sigmoidabTransform(1/(2*(lengthscale^2)), 'xtoa', inversewidth_range);  
  end;
  
  % initialization of POL2 effect variance
  pars(pol2effectvar_index)=sigmoidabTransform(0.7*var(dataVals1), 'xtoa', pol2effectvar_range);
  % initialization of POL2 noise variance
  pars(pol2noisevar_index)=sigmoidabTransform(0.1*var(dataVals1), 'xtoa', pol2noisevar_range);
end;

if (initializationtype>2),
%  if isempty(lengthscale),
%    % initialization of inverse squared width
%    pars(inversewidth_index)=log(1/(2*(15^2)));
%  else
%    pars(inversewidth_index)=log(1/(2*(lengthscale^2)));  
%  end;

randval1=rand()
randval2=rand()
randval3=rand()
randval4=rand()
randval5=rand()
randval6=rand()
randval7=rand()
randval8=rand()
randval9=rand()
randval10=rand()

  tempscale=(max(timevector)-min(timevector))*randval1;
  pars(inversewidth_index)=sigmoidabTransform(1/(2*(tempscale^2)), 'xtoa', inversewidth_range);  


  % initialization of POL2 effect variance
  pars(pol2effectvar_index)=sigmoidabTransform(randval2*var(dataVals1), 'xtoa', pol2effectvar_range);
  % initialization of POL2 noise variance
  pars(pol2noisevar_index)=sigmoidabTransform(randval3*var(dataVals1), 'xtoa', pol2noisevar_range);
end;

%pars
pol2model = gpdisimExpandParam(pol2model,pars);
