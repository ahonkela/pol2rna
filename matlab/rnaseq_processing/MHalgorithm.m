function[THETA mu1 ACC]=MHalgorithm(a,t1_mcmc,t2_mcmc,t3_mcmc,no_iter)

% a can be at most 78! (we have 78 gene groups)
   
m_combined=[mean(log(t1_mcmc),2) mean(log(t2_mcmc),2) mean(log(t3_mcmc),2)];
m=mean(m_combined,2);
[Y0 I0]=sort(m);


ii=[110;120;130;140;150;160;170;180;190;200];

if a<78
    A=(a-1)*500+1;
    B=a*500;
else 
    A=38501;
    B=39071;
end

ind1=I0(A:B); % first group of transcripts; indices
mu1=Y0(A:B); % ; combined mean

mu0=log(t1_mcmc(ind1,:))'; 
mu5=log(t2_mcmc(ind1,:))';
mu10=log(t3_mcmc(ind1,:))';

alpha_start=1;
beta_start=1;
theta_start=[alpha_start; beta_start];

ALPHA=zeros(no_iter,500);
BETA=zeros(no_iter,500);
ACC=zeros(500,1); 
ALPHA_thinned=zeros(10,500);
BETA_thinned=zeros(10,500);
meanALPHA=zeros(500,1);
meanBETA=zeros(500,1);
meanParams=zeros(2,500);

for k=1:500
    
	y=[mu0(k,:)' mu5(k,:)' mu10(k,:)']; %500*3
	theta_hat=param_map_estimate(@(th) negloglike(th,y), theta_start); %MAP estimate
	H=compute_hessian(theta_hat,y);
	Sigma_propose=((2.38^2)/2)*inv(-H);

	theta=[];
	acc=0;
	llik=-negloglike(theta_hat,y);
	ll=llik;
	theta_t=theta_hat;
	negloglike_t=-ll;
	R_propose=jitChol(Sigma_propose);

	for t=1:no_iter
		theta_prop=abs(theta_t+R_propose'*randn(length(theta_t),1));
		llik=-negloglike(theta_prop,y);
    		ll_prop=llik;
    		negloglike_prop=-ll_prop;
    		acp_ratio=exp(-negloglike_prop+negloglike_t);
		if acp_ratio>1
			% accept
			theta_t = theta_prop;
			negloglike_t = negloglike_prop;
			if t>100
				acc=acc+1;
			end
		else
		r=rand;
			if r<acp_ratio
				%accept
				theta_t = theta_prop;
				negloglike_t = negloglike_prop;
				if t>100
					acc=acc+1;
				end
			end
		end
		theta= [theta theta_t];
	end
	ALPHA(:,k)=theta(1,:);
	BETA(:,k)=theta(2,:);
	acceptance_rate=acc/(T/2);
	ACC(k)=acceptance_rate;

	ALPHA_thinned(:,k)=ALPHA(ii,k);
	BETA_thinned(:,k)=BETA(ii,k);
	meanALPHA(k)=mean(ALPHA_thinned(:,k));
	meanBETA(k)=mean(BETA_thinned(:,k));
	meanParams(:,k)=[meanALPHA(k); meanBETA(k)];

end

THETA=mean(meanParams,2);

end

