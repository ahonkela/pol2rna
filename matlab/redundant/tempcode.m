initializationtype=1;
[model1,temptransforminfo1]=createSimDisim(temptimes,tempvals1,tempvals2,lengthscale,initializationtype);
[pars1,nams1]=gpdisimExtractParam(model1);
parst1=transformParametersWithSettings(pars1,temptransforminfo1,'atox');

parst2=parst1;
parst2(4)=1e-1;
parst2(5)=1e-2;
parst2(7)=80;
parst2(8)=1e-3;
parst2(9)=1e-3;
pars2=transformParametersWithSettings(parst2,temptransforminfo1,'xtoa');
model2=gpdisimExpandParam(model1,pars2);
temptimes=[0:20:1280]';
[tempvals1,tempvals2,temptimes]=generateDataFromModel(model2,temptimes);
clf;plot(sqrt(temptimes),tempvals1,'r-');hold on;plot(sqrt(temptimes),tempvals2,'g-');

maxiters=100; ninits=10; lengthscale=5;
[model3,temptransforminfo3,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createGeneGPModels(temptimes,tempvals1,tempvals2,lengthscale,maxiters,ninits);
[pars3,nams3]=gpdisimExtractParam(model3);
parst3=transformParametersWithSettings(pars3,temptransforminfo3,'atox');

plotpredictions(model3,[0:5:1280]',2,1,1,'exampletitle');


maxiters=100; ninits=10; lengthscale=5;
[jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createGeneGPModels(temptimes,tempvals1,tempvals2,lengthscale,maxiters,ninits);
model1=jointmodelb;
temptransforminfo1=jointtransforminfo;
[pars1,nams1]=gpdisimExtractParam(model1);
parst1=transformParametersWithSettings(pars1,temptransforminfo1,'atox');
parst2=parst1;
parst2(11)=0.695;
parst2(10)=0.0077;
parst2(4)=0.0128;
parst2(5)=0.0103^2;
parst2(7)=120;
parst2(8)=1e-3;
parst2(9)=1e-3;
pars2=transformParametersWithSettings(parst2,temptransforminfo1,'xtoa');
model2=gpdisimExpandParam(model1,pars2);
predtimes=[0:5:1280]';
[pol2preds,pol2vars,rnapreds,rnavars]=plotpredictions(model2,predtimes,2,0.01,0,'exampletitle');

model3=gpdisimOptimise(model2);
[pars3,nams3]=gpdisimExtractParam(model3);
parst3=transformParametersWithSettings(pars3,temptransforminfo1,'atox');



rnaest=zeros(size(pol2preds));
for k=1:length(rnafrompol),
  delayedt=max([predtimes(k)-delay 0]);
  rnaest(k)=rna0*exp(-D*predtimes(k))...
	    +B/D*(1-exp(-D*predtimes(k)))...
	    +polmean*S/D*(1-exp(-D*delayedt));
  if k==1,
    tempint==0;
  else
    timesteps=predtimes(2:k)-predtimes(1:k-1);

    delayedtimes=predtimes(1:k-1)-delay;
    delayedindices=nan*delayedtimes;
    delayedpreds=nan*delayedtimes;
    for l=1:delayedtimes,
      if delayedtimes(l) >= 0,
      tempind=max(find(predtimes<=delayedtimes(l)));
      
      delayedindices=max(find(predtimes<=
      % use linear interpolation
    end;
    
      
    
    tempint=sum(timesteps.*exp(D*predtimes(1:k-1)).*pol2preds
  rnaest(k)=rnaest(k)+S*exp(-D*t)
	    
end;



model=results_jointmodels{k};
trinfo=results_jointtransforminfos{k};
ll=results_loglikelihoods(k,3);
ensemblid=results_ensemblids(k);
[pars,nams]=gpdisimExtractParam(model)
parst=transformParametersWithSettings(pars,trinfo,'atox');
mytitle=sprintf('ENSG %d (Pol2+RNA), LL %f, m(0)=%f, B=%f, S=%f, \\Delta=%f', ...
                ensemblid,ll,parst(11),parst(10),sqrt(parst(5)),parst(4));
plotpredictions(results_jointmodels{1},[0:5:1280]',2,1,3,mytitle);
myfilename=sprintf('simdisim_polrna_ensg%d.eps',ensemblid);
print(myfilename,'-depsc');
