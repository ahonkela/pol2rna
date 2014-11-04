function[m_groups biol_var]=estimate_biolvar(t1_mcmc,t2_mcmc,t3_mcmc)

alpha=zeros(78,1);
beta=zeros(78,1);
v=zeros(78,1);
m_groups=zeros(78,1);
stdev=zeros(78,1);
acc=zeros(78,1);

no_iter=200;

for i=1:78
	[THETA mu1 ACC]=MHalgorithm(i,t1_mcmc,t2_mcmc,t3_mcmc,no_iter);
	alpha(i)=THETA(1);
	beta(i)=THETA(2);
	v(i)=beta(i)/alpha(i);
	m_groups(i)=mean(mu1);
	acc(i)=mean(ACC);
end

biol_var = smooth(m_groups,v,0.1,'loess');

end
