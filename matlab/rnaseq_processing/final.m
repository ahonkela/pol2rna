% Inputs:
% - mu: 39071-by-10 matrix which contains the mean log-rpkm gene expression levels
% - V: 39071-by-10 matrix which contains the technical variances of the log-rpkm gene expression levels
% - t1_mcmc: 39071-by-500 matrix which contains all BitSeq MCMC samples for time point 0 (non-logged).
% - t2_mcmc: 39071-by-500 matrix which contains all BitSeq MCMC samples for time point 5 (non-logged).
% - t3_mcmc: 39071-by-500 matrix which contains all BitSeq MCMC samples for time point 10 (non-logged).

function[]=final(mu,V,t1_mcmc,t2_mcmc,t3_mcmc)

[m_groups biol_var]=estimate_biolvar(t1_mcmc,t2_mcmc,t3_mcmc);
GPanalysis(mu,V,m_groups,biol_var);

end
