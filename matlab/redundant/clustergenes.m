addpath /share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/pol2rnaseq/matlab
load /share/synergy/analyses/all_profiles_2012-11-21.mat
cd /share/mi/workspace/jtpelto/synergy/synergy_data/analyses/clustering_2012_11_26

% Normalize profiles to be zero mean and unit variance in both POL2 and mRNA 
munormalized=mu;
for k=1:size(munormalized,1),
  temppol = mu(k,[1:101]);
  temprna = mu(k,[102:202]);
  munormalized(k,:) = [(temppol-mean(temppol))/sqrt(var(temppol)) (temprna-mean(temprna))/sqrt(var(temprna))];
end;

% Divide data into a training and test set, reduce feature dimensionality
randvector = rand(size(munormalized,1),1);
trainingsetI = find(randvector<0.5);
testsetI = find(randvector >= 0.5);
featuresetI = [[1:5:101] [102:5:202]];
mutrain = munormalized(trainingsetI,featuresetI);
mutest = munormalized(testsetI,featuresetI);


% For each potential cluster number, train the GMM, then compute log-likelihood of held-out data
niter = 100;

maxclust=50;
loglikelihoods=zeros(maxclust,2);
for nclust=1:maxclust,
  [clustmemberships,clustmeans,clustcovs,clustprior,totalloglikelihoodtrain] = gaussmixt(mutrain, ones(size(mutrain,1),1), nclust, niter, [],[],[]);

  [clustmemberships2,clustmeans2,clustcovs2,clustprior2,totalloglikelihoodtest] = gaussmixt(mutest, ones(size(mutest,1),1), nclust, 0, clustmeans,clustcovs,clustprior);

  loglikelihoods(nclust,:) = [totalloglikelihoodtrain totalloglikelihoodtest];
  fprintf(1,'Loglikelihoods at %d clusters: %e (train) %e (test)\n', nclust, totalloglikelihoodtrain, totalloglikelihoodtest);
end;

% Log-likelihood plot
clf; plot(loglikelihoods(:,1),'r-'); hold on; plot(loglikelihoods(:,2),'g-');
xlabel('Number of clusters'); ylabel('RED=train-likelihood,GREEN=test-likelihood');

% BIC plot
% clf; k0=[1:172]';k=k0.*(42 + 0*(42*41)/2+42 + (k0-1)./k0); plot(-2*loglikelihoods(:,1)+k*log(494),'r-'); hold on; plot(-2*loglikelihoods(:,2)+k*log(522),'g-');



bestnclust=20;
nrepeats=50,
loglikelihoods2=zeros(nrepeats,1);
clusterings=cell(nrepeats,4);
for irepeat=1:nrepeats,
  [clustmemberships,clustmeans,clustcovs,clustprior,totalloglikelihoodtrain] = gaussmixt(munormalized(:,featuresetI), ones(size(munormalized,1),1), bestnclust, niter, [],[],[]);
  loglikelihoods2(irepeat)=totalloglikelihoodtrain;
  clusterings{irepeat,1} = clustmemberships;
  clusterings{irepeat,2} = clustmeans;
  clusterings{irepeat,3} = clustcovs;
  clusterings{irepeat,4} = clustprior;
end;

Ibest = find(loglikelihoods2==max(loglikelihoods2));
clustmemberships=clusterings{Ibest,1};
clustmeans=clusterings{Ibest,2};
clustcovs=clusterings{Ibest,3};
clustprior=clusterings{Ibest,4};


hardmemberships = zeros(size(munormalized,1),1);
for k=1:size(munormalized,1),
    I = find(clustmemberships(k,:)==max(clustmemberships(k,:)));
    hardmemberships(k) = I;
end;

hardclustersizes=zeros(bestnclust,1);
hardclustmeans=zeros(bestnclust,size(munormalized,2));
hardclustvars=zeros(bestnclust,size(munormalized,2));
for k=1:bestnclust,
    Ik = find(hardmemberships==k);
    hardclustersizes(k)=length(Ik);
    hardclustmeans(k,:) = mean(munormalized(Ik,:));
    hardclustvars(k,:) = var(munormalized(Ik,:));
end;

featuresetIPol2=[1:5:101];
featuresetIRNA=[102:5:202];
timeindices=[1:21]/22*1280;
for k=1:bestnclust,
%    clf; errorbar([[1:5:101] [102:5:202]],hardclustmeans(k,featuresetI),sqrt(hardclustvars(k,featuresetI)));
    clf; errorbar(timeindices,hardclustmeans(k,featuresetIPol2),sqrt(hardclustvars(k,featuresetIPol2)),'r-');
    hold on; errorbar(timeindices+2,hardclustmeans(k,featuresetIRNA),sqrt(hardclustvars(k,featuresetIRNA)),'g-');

    title(sprintf('Cluster %d (%d members)\n',k,hardclustersizes(k)));
    print(sprintf('cluster%d_corrected.eps',k),'-depsc');
end;

save hardmemberships_corrected_2012_11_26.txt hardmemberships -ASCII




genenames=read_stringfile('~/tigretemp/hmcplots_2012_11_21/hmc_delay_prctiles_2012-11-21.txt',[' ' 10 13]);
f = fopen('clusterindices_corrected_2012_11_26.txt','w');
for k=1:size(hardmemberships,1),
  fprintf(f,'%s %d\n', genenames{k+1}{1}, hardmemberships(k));
end;
fclose(f);
