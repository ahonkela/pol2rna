function [rnapred,rnapred2]=manualRNAPred(pol2pred, times, startmean,decay,variance, ...
			    basal);

sensitivity=sqrt(variance);

rnapred=zeros(size(pol2pred));
rnapred(1)=startmean;
for k=2:length(times),
  timestep=times(k)-times(k-1);
  rnapred(k)=rnapred(k-1)+timestep*basal...
      +timestep*(sensitivity*pol2pred(k-1)-decay*rnapred(k-1));
end;




rnapred2=zeros(size(pol2pred));
rnapred2(1)=startmean;
for k=2:length(times),
  timesteps=times(2:k)-times(1:k-1);
  tempsum=sum(timesteps.*exp(decay*times(1:k-1)).*pol2pred(1:k-1));
  rnapred2(k)=startmean*exp(-decay*times(k))+basal/decay*(1-exp(-decay*times(k)))...
      +sensitivity*exp(-decay*times(k))*tempsum;
end;

