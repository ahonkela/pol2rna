% compute kernel between asynchronously observed POL2 values
% (actually, smoothed POL2 values where the smoothing comes from
% increasing asynchrony of several cells that are averaged to get
% the observation).

function kernelvalue=kernel_async_p(t,tb,A,L,sigma2);

kernelvalue=A*L/sqrt(L^2 + 2*sigma2*tb + 2*sigma2*t)*exp(-(t-tb)^2/(L^2 + 2*sigma2*tb + 2*sigma2*t));

