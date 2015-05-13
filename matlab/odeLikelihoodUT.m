function ll = odeLikelihoodUT(p, m, mvar, t_obs, params, doplot, includeprior),

if nargin < 6,
  doplot = 0;
end

if nargin < 7,
  includeprior = 0;
end

[tout, yout] = ode45(@(t, y) delayode(t, y, params, p, t_obs), [t_obs], params(5));
ll = gaussianLogLikelihood(yout, m, mvar);
if includeprior,
  ll = ll + odeLogPrior(params);
end

if doplot,
  [tout, yout] = ode45(@(t, y) delayode(t, y, params, p, t_obs), [min(t_obs), max(t_obs)], params(5));
  subplot(2, 1, 1);
  plot(t_obs-min(t_obs), p, 'rx');
  v = axis;
  axis([0, 1300, 0, v(4)])
  subplot(2, 1, 2);
  plot(tout-min(t_obs), yout);
  hold on
  plot(t_obs-min(t_obs), m, 'gx');
  hold off
  v = axis;
  axis([0, 1300, 0, v(4)])
end


function grad = delayode(t, m, params, p, t_obs),

beta0 = params(1);
beta2 = params(2);
beta = sqrt(beta2);
%alpha = log(2)/params(3);
alpha = params(3);
Delta = params(4);

grad = beta0 + beta*joinTheDots(p, t_obs, t-Delta) - alpha * m;


function val = joinTheDots(p, t_obs, t),

if t < t_obs(1),
  val = p(1);
else
  val = interp1(t_obs, p, t, 'linear');
end


function ll = gaussianLogLikelihood(x, y, vars),

ll = sum(-.5*numel(x)*log(2*pi) - .5*sum(log(vars)) ...
	 -.5*sum((x-y).^2 ./ vars));


function ll = odeLogPrior(params),

beta0 = params(1);
beta2 = params(2);
beta = sqrt(beta2);
%alpha = log(2)/params(3);
alpha = params(3);
Delta = params(4);
m0 = params(5);

ll = lnlogitnormalpdf(beta0, 0, 2, 0, 1) + ...
     lnlogitnormalpdf(beta2, 0, 2, 1e-6, 1) + ...
     lnlogitnormalpdf(alpha, 0, 2, 1e-6, log(2)) + ...
     lnlogitnormalpdf(Delta, -2, 2, 0, 300) + ...
     lnlogitnormalpdf(m0, 0, 2, 0, 10);
