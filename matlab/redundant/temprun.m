for k=kstart:kend,
  gene_index=interestinggenes(k);  
  randn('seed',gene_index);
  rand('seed',gene_index+1234567);
  
  dataVals1=pol_summaryseries(gene_index,:)';
  dataVals2=rna_summaryseries(gene_index,:)';
  timevector=measurementtimes'; 
  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;
  
  maxiters=100; ninits=8; lengthscale=0.25;
    
  % set parameter ranges
  inversewidth_range=1./([1280 0.1].^2);
  pol2effectvar_range=[0.0002 20];%*var(tempvals1);
  if pol2effectvar_range(2)==pol2effectvar_range(1),
    pol2effectvar_range(2)=pol2effectvar_range(1)+1;
  end;  
  pol2noisevar_range=[0.005 5];%*var(tempvals1);
  if pol2noisevar_range(2)==pol2noisevar_range(1),
    pol2noisevar_range(2)=pol2noisevar_range(1)+1;
  end;  
  rnadecay_range=[0.00001 1.5];
  rnaeffectvar_range=[0 20];
  rnadelay_range=[0 5];
  rnanoisevar_range=[0.005 5];%*var(tempvals2);
  if rnanoisevar_range(2)==rnanoisevar_range(1),
    rnanoisevar_range(2)=rnanoisevar_range(1)+1;
  end;
  rnabasal_range=[0 10000];
  rnastartmean_range=[0 10000];
  pol2mean_range=[0 10000];
  parameterranges=[inversewidth_range;pol2effectvar_range;pol2noisevar_range;rnadecay_range;rnaeffectvar_range;rnadelay_range;rnanoisevar_range;rnabasal_range;rnastartmean_range;pol2mean_range];  
      
  use_fixedrnavariance=0;
  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createNdGeneGPModels(temptimes,tempvals1,tempvals2,lengthscale,maxiters,ninits,parameterranges,use_fixedrnavariance);

  fprintf(1, 'Gene %d (name %s, nanoid %s): optimized loglik: %f\n',...
	  k,genenames{interestinggenes(k)},genenanoids{interestinggenes(k)},joint_ll);

  results_geneindices(k)=gene_index;
  results_genenames{k}=genenames{gene_index};
  results_loglikelihoods(k,:)=[naive_ll rbf_ll joint_ll];
  results_jointmodels{k}=jointmodelb;
  results_jointtransforminfos{k}=jointtransforminfo;
  if mod(k,5)==0,
    tempfilename=sprintf('fittingresults_rnaurna_temp_%d.mat',startpercent);
    save(tempfilename,'results_geneindices','results_genenames',...
       'results_loglikelihoods','results_jointmodels','results_jointtransforminfos','-mat');
  end;
end;

tempfilename=sprintf('fittingresults_rnaurna_temp_%d.mat',startpercent);
save(tempfilename,'results_geneindices','results_genenames',...
     'results_loglikelihoods','results_jointmodels','results_jointtransforminfos','-mat');
