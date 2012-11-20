load all_profiles_2012-10-08b.mat

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

maxclust=200;
loglikelihoods=zeros(maxclust,2);
for nclust=68:maxclust,
  [clustmemberships,clustmeans,clustcovs,clustprior,totalloglikelihoodtrain] = gaussmixt(mutrain, ones(size(mutrain,1),1), nclust, niter, [],[],[]);

  [clustmemberships2,clustmeans2,clustcovs2,clustprior2,totalloglikelihoodtest] = gaussmixt(mutest, ones(size(mutest,1),1), nclust, 0, clustmeans,clustcovs,clustprior);

  loglikelihoods(nclust,:) = [totalloglikelihoodtrain totalloglikelihoodtest];
  fprintf(1,'Loglikelihoods at %d clusters: %e (train) %e (test)\n', nclust, totalloglikelihoodtrain, totalloglikelihoodtest);
end;

% Log-likelihood plot
clf; plot(loglikelihoods(:,1),'r-'); hold on; plot(loglikelihoods(:,2),'g-');

% BIC plot
% clf; k0=[1:172]';k=k0.*(42 + 0*(42*41)/2+42 + (k0-1)./k0); plot(-2*loglikelihoods(:,1)+k*log(494),'r-'); hold on; plot(-2*loglikelihoods(:,2)+k*log(522),'g-');

bestnclust=20;
  [clustmemberships,clustmeans,clustcovs,clustprior,totalloglikelihoodtrain] = gaussmixt(munormalized(:,featuresetI), ones(size(munormalized,1),1), bestnclust, niter, [],[],[]);



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

for k=1:bestnclust,
    clf; errorbar([[1:5:101] [102:5:202]],hardclustmeans(k,featuresetI),sqrt(hardclustvars(k,featuresetI)));
    title(sprintf('Cluster %d (%d members)\n',k,hardclustersizes(k)));
    print(sprintf('cluster%d_corrected.eps',k),'-depsc');
end;

save hardmemberships_corrected.txt hardmemberships -ASCII




genenames=read_stringfile('~/tigretemp/hmcplots_2012_10_08/hmc_delay_prctiles_2012-10-08b.txt',[' ' 10 13]);
f = fopen('clusterindices_corrected.txt','w');
for k=1:size(hardmemberships,1),
  fprintf(f,'%s %d\n', genenames{k+1}{1}, hardmemberships(k));
end;
fclose(f);
