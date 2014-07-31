function [clustmemberships,clustmeans,clustcovs,clustprior,totalloglikelihood] = gaussmixt(x, xweights, nclust, niter, initmeans, initcovs, initprior);
%---------------------------------------------------------
% Uses the Expectation-Maximization (EM) algorithm to fit a Gaussian mixture model to the data.
% Implementation Jaakko Peltonen, October 2012
%
% Inputs:
% x : data matrix, ndata*ndim
% xweights : weights for data points, they do not currently affect the code.
% nclust : desired number of clusters
% niter : number of EM iterations to run
% initmeans : optional parameter, nclust*ndim matrix, mean vector of each cluster
% initcovs : optional parameter, nclust*ndim*ndim matrix, covariance matrix of each cluster
% initprior : optional parameter, nclust*1 vector, prior probability of each cluster, must sum to 1
%
% Outputs:
% clustmemberships : ndata*nclust matrix, membership of each data point in each cluster (sums to 1)
% clustmeans : nclust*ndim matrix, mean vector of each cluster
% clustcovs : nclust*ndim*ndim matrix, covariance matrix of each cluster
% clustprior : nclust*1 vector, prior probability of each cluster
% totalloglikelihood : scalar, log-likelihood of the fitted mixture model (log-probability of the data)
%---------------------------------------------------------

drawclusters=0;

ndata=size(x,1);
ndim=size(x,2);

% Initialize cluster means at given initialization or at random data points
if ~isempty(initmeans),
  clustmeans = initmeans;
else
  clustmeans=zeros(nclust,ndim);
  for k=1:nclust,
      l=ceil(ndata*rand);
      clustmeans(k,:)=x(l,:);
  end;
end;

% Initialize covariance matrices at given initialization or to small diagonal values
if ~isempty(initcovs),
  clustcovs = initcovs;
else
  clustcovs=zeros(nclust,ndim,ndim);
  smallcov=(1e-6)^(1/ndim);
  for k=1:nclust,
      clustcovs(k,:,:)=eye(ndim)*smallcov;
  end;
end;

% Initialize cluster prior probabilities to given initialization or to a uniform distribution
if ~isempty(initprior),
  clustprior = initprior;
else
  clustprior=ones(nclust,1)/nclust;
end;

clustdist=zeros(ndata,nclust);
clustinvcovs=zeros(nclust,ndim,ndim);
clustdets=zeros(nclust,1);
clustmemberships=zeros(ndata,nclust);

for i=1:niter+1,
    % Compute inverse covariances and determinants of clusters
    for k=1:nclust,
      clustinvcovs(k,:,:)=inv(squeeze(clustcovs(k,:,:)));
      clustdets(k)=det(squeeze(clustcovs(k,:,:)));
    end;
    clustcovdiags=zeros(nclust,ndim);
    for k=1:nclust,
      clustcovdiags(k,:)=diag(squeeze(clustcovs(k,:,:)));
    end;
    % log(clustcovdiags)
    % clustdets
    % pause

    %------------------------------
    % E-STEP
    %------------------------------
    % Compute Mahalanobis distances to clusters, 
    % assign unnormalized cluster memberships (log-likelihoods of data)
    for k=1:nclust,
        tempdiff = (x-repmat(clustmeans(k,:),[ndata 1]));
	tempdiff2 = tempdiff*squeeze(clustinvcovs(k,:,:));
        clustdist(:,k) = sum(tempdiff.*tempdiff2,2);
	clustmemberships(:,k) = -clustdist(:,k)/2 -(ndim/2)*log(2*pi) - (1/2)*log(clustdets(k)) +log(clustprior(k));
    end;
    %clustdist
    %imagesc(clustdist)
    %drawnow
    %pause
    % clustmemberships
    % pause

    % Assign points to clusters by normalizing the assignments; 
    % compute data log-likelihood at the same time
    clustmaxmemberships = max(clustmemberships,[],2);
    clustmemberships = exp(clustmemberships - repmat(clustmaxmemberships,[1 nclust]));

    dataloglikelihoods = log(sum(clustmemberships,2)) + clustmaxmemberships;
    totalloglikelihood = sum(dataloglikelihoods);
    fprintf(1,'iteration %d, totalloglikelihood %e\n',i, totalloglikelihood);
    % pause

    clustmemberships = clustmemberships./repmat(sum(clustmemberships,2),[1 nclust]);
    % sum(sum(clustmemberships))
    if (drawclusters)
      imagesc(clustmemberships);
      drawnow;
    end;
    % pause

    if (i < niter),
      %------------------------------
      % M-STEP
      %------------------------------
      % Update cluster prior
      clustprior = sum(clustmemberships,1)'/ndata;
    
      for k=1:nclust,
        if clustprior(k)>0,
          % Update cluster means
          clustmeans(k,:) = sum(x.*repmat(clustmemberships(:,k),[1 ndim]),1)/sum(clustmemberships(:,k));
  
          % Update cluster covariance matrices
          tempcov = zeros(ndim,ndim);
          for l=1:ndata,
            tempcov = tempcov + clustmemberships(l,k)*((x(l,:)-clustmeans(k,:))')*(x(l,:)-clustmeans(k,:));
          end;
          tempcov = tempcov/sum(clustmemberships(:,k));

          % regularization
          smallcov=(1e-6)^(1/ndim);
          tempcov=tempcov + eye(ndim)*smallcov;

          clustcovs(k,:,:)=tempcov;
        else
          % Reinitialize cluster mean at a random data point
          l=ceil(ndata*rand);
          clustmeans(k,:)=x(l,:);
          % Reinitialize covariance matrix to small diagonal values
          smallcov=(1e-6)^(1/ndim);
          clustcovs(k,:,:)=eye(ndim)*smallcov;
          % Reinitialize cluster prior probabilities to a value corresponding to a uniform distribution
          clustprior(k)=1/nclust;
        end;

      end;

      % renormalize cluster prior in case any clusters were reinitialized
      clustprior = clustprior/sum(clustprior);
    end;
end;

