function [clustind,clustmeans,clustcovs] = kmeans(x, xweights, nclust, niter);

ndata=size(x,1);
ndim=size(x,2);

% Initialize cluster means at random data points
clustmeans=zeros(nclust,ndim);
for k=1:nclust,
    l=ceil(ndata*rand);
    clustmeans(k,:)=x(l,:);
end;


clustdist=zeros(ndata,nclust);
for i=1:niter,
    i
    % Compute distances to clusters
    for k=1:nclust,
        clustdist(:,k) = sum((x-repmat(clustmeans(k,:),[ndata 1])).^2,2);
    end;
    % Assign points to closest cluster
    [mindist,clustind] = min(clustdist');
    sum(mindist)
    clustind=clustind';
    size(clustind)
    size(xweights)
    
    % Update cluster means, taking data weights into account
    for k=1:nclust,
        % Update non-empty clusters, reinitialize empty ones at
        % random data points.
        if sum((clustind==k).*xweights)==0,
            l=ceil(ndata*rand);
            clustmeans(k,:)=x(l,:);
        else
            clustmeans(k,:)=sum(repmat((clustind==k).*xweights,[1 ndim]).*x,1)/sum((clustind==k).*xweights);
        end;
    end;
end;

% calculate covariances around clusters
clustcovs=cell(nclust,1);
for k=1:nclust,
    if sum((clustind==k).*xweights)==0,
        clustcovs{k} = eye(ndim);
    else
        tempcov=zeros(ndim,ndim);
        for l=1:ndim,
            for m=1:ndim,
                tempcov(l,m) = sum((x(:,l)-repmat(clustmeans(k,l),[ndata 1])).*(x(:,m)-repmat(clustmeans(k,m),[ndata 1])).*xweights.*(clustind==k))/sum((clustind==k).*xweights);
            end;
        end;
        clustcovs{k} = tempcov;
    end;
end;
