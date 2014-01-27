function [erfdiff,logerfdiff,erfdsign]=myerfdiff(x1,x2)
% computes erf(x1)-erf(x2) by slow numerical integration


nintervals=10000;
interval=(x1-x2)/nintervals;

scaling=max([-x1*x1 -x2*x2]);


erfdiff=0;
for k=0:(nintervals-1),
  tempx=x2+k*interval;
  erfdiff=erfdiff+exp(-tempx*tempx-scaling);
end;
erfdiff
erfdsign=sign(erfdiff)*sign(interval)

logerfdiff=log(abs(erfdiff)) + log(abs(interval)) + scaling + log(2) -0.5*log(pi);
erfdiff=erfdsign*exp(logerfdiff);


