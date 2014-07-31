function [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createNdGeneGPModels_celltimes_newdata(timevector,dataVals,lengthscale,maxiters,ninits,parameterranges,use_fixedrnavariance);


DO_PLOT = 0;


% Naive log-likelihood of POL2
dataVals1=dataVals{1};
if isempty(dataVals1)==0,
  datatemp=dataVals1-min(dataVals1);
  if mean(datatemp.^2)>0,
    datatemp=datatemp/sqrt(mean(datatemp.^2));
    dim = size(datatemp, 1);
    tempK=eye(dim)*var(datatemp);
    tempm=datatemp-mean(datatemp);
    ll1 = -(dim/2)*log(2*pi) - (1/2)*log(det(tempK)) - (1/2)*tempm'*inv(tempK)*tempm;
  else
    ll1=inf;
  end;
else
  ll1 = 0;
end;

% Naive log-likelihood of RNA
if (length(dataVals)>1) && (isempty(dataVals{2})==0),
  dataVals2=dataVals{2};
  if mean(dataVals2.^2)>0,
    datatemp=dataVals2/sqrt(mean(dataVals2.^2));
    dim = size(datatemp, 1);
    tempK=eye(dim)*var(datatemp);
    tempm=datatemp-mean(datatemp);
    ll2 = -(dim/2)*log(2*pi) - (1/2)*log(det(tempK)) - (1/2)*tempm'*inv(tempK)*tempm;
  else
    ll2=inf;
  end;
else
  ll2 = 0;
end;

% Naive log-likelihood of independently modeled POL2 and RNA
naive_ll=ll1+ll2
%pause

create_rbfmodels=0;

if create_rbfmodels==1,
  %--------------
  % Independent ARBF models for both POL2 and RNA
  %--------------
  pol2_ll=-inf;
  for initializationtype=1:ninits,
    pol2model = createRBFModel(timevector,dataVals1,lengthscale,initializationtype);
    [pars,nams]=gpdisimExtractParam(pol2model);
    %ninits=5;
    temppars=pars;
    temppars
    tempmodel=gpdisimExpandParam(pol2model,temppars);
    fprintf(1,'POL2 rbf model numGenes %d\n',tempmodel.numGenes);
    tempmodelb = gpdisimOptimise(tempmodel,1,maxiters);
    templl=gpdisimLogLikelihood(tempmodelb);
    % pause
    if templl>pol2_ll,
      pol2modelb=tempmodelb;
      pol2_ll=templl;
    end;  
  end;
  
  
  rna_ll=-inf;
  for initializationtype=1:ninits,
    rnamodel = createRBFModel(timevector,dataVals2,lengthscale,initializationtype);
    [pars,nams]=gpdisimExtractParam(rnamodel);
    %ninits=5;
    temppars=pars;
    tempmodel=gpdisimExpandParam(rnamodel,temppars);
    tempmodelb = gpdisimOptimise(tempmodel,1,maxiters);
    templl=gpdisimLogLikelihood(tempmodelb);
    if templl>rna_ll,
      rnamodelb=tempmodelb;
      rna_ll=templl;
    end;  
  end;

  rbf_ll=pol2_ll+rna_ll
  %pause
else
  pol2modelb=nan;
  rnamodelb=nan;
  rbf_ll=nan;
end;







joint_ll=-inf;
bestinit=-1;
addPriors = 0;
for initializationtype=1:ninits,
  %try
    fprintf(1,'Trying to create initialization %d for SIM-DISIM\n',initializationtype);
    [jointmodel,temptransforminfo]=createNdSimDisim_celltimes_newdata(timevector,dataVals,lengthscale,initializationtype,parameterranges,use_fixedrnavariance,addPriors);
    fprintf(1,'Trying to create initialization %d for SIM-DISIM, step2\n',initializationtype);
    init_ll=gpnddisimLogLikelihood(jointmodel);
  
    fprintf(1,'Trying to create initialization %d for SIM-DISIM, initll=%f, step3\n',initializationtype,init_ll);
    [pars,nams]=gpnddisimExtractParam(jointmodel)


if DO_PLOT,  
    ensemblid=186628; timeshift=300; predicttimes=timeshift + 1280*(([0:256]'/256).^2);
    plotpredictions(jointmodel,predicttimes,2,1,0,1,1,0,[],timeshift);  
    lltemp=gpnddisimLogLikelihood(jointmodel);
    plottitle=makeplottitle(jointmodel,lltemp,ensemblid);
    plottitle=[sprintf('INIT%d (iWidth %e,var %e) ',initializationtype,jointmodel.kern.comp{1}.comp{1}.inverseWidth,jointmodel.kern.comp{1}.comp{1}.variance) plottitle];
    axes; axis([0 1 0 1]); h=text(0,1.07,plottitle); 
    set(h,'fontname','Helvetica'); set(h,'fontsize',20); axis off;
    myfilename=sprintf('ENSG%011d_pol2_init%d.eps',ensemblid,initializationtype)
    print(myfilename,'-depsc','-S700,600');
    %pause
end;
        
    fprintf(1,'Trying to optimize from initialization %d for SIM-DISIM\n',initializationtype);
  %  try
      tempmodelb = gpnddisimOptimise(jointmodel,1,maxiters);
  %  catch
  %    tempmodelb=jointmodel;
  %  end;


        
    fprintf(1,'Trying to evaluate model trained from initialization %d for SIM-DISIM\n',initializationtype);
    templl=gpnddisimLogLikelihood(tempmodelb);
    fprintf(1,'Initialization %d: init_ll %f, optimized %f\n', initializationtype, init_ll, templl);
    %pause      
    fprintf(1,'Done processing initialization %d for SIM-DISIM, best init=%d\n',initializationtype,bestinit);
 % catch
 %   templl=-inf;
 % end;

if DO_PLOT,
    ensemblid=186628; timeshift=300; predicttimes=timeshift + 1280*(([0:256]'/256).^2);
    plotpredictions(tempmodelb,predicttimes,2,1,0,1,1,0,[],timeshift);  
    lltemp=gpnddisimLogLikelihood(tempmodelb);
    plottitle=makeplottitle(tempmodelb,lltemp,ensemblid);
    plottitle=[sprintf('TRAINED%d (iWidth %e,var %e) ',initializationtype,tempmodelb.kern.comp{1}.comp{1}.inverseWidth,tempmodelb.kern.comp{1}.comp{1}.variance) plottitle];
    axes; axis([0 1 0 1]); h=text(0,1.07,plottitle); 
    set(h,'fontname','Helvetica'); set(h,'fontsize',20); axis off;
    myfilename=sprintf('ENSG%011d_pol2_trained%d.eps',ensemblid,initializationtype)
    print(myfilename,'-depsc','-S700,600');
    %pause
end;

    
  if templl>joint_ll,
    jointmodelb=tempmodelb;
    jointtransforminfo=temptransforminfo;
    joint_ll=templl;
    bestinit=initializationtype;
  end;  
end;




