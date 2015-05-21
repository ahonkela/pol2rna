function ll = odeLikelihoodSpline(p, m, mvar, t_obs, params, doplot, includeprior),

if nargin < 6,
  doplot = 0;
end

if nargin < 7,
  includeprior = 0;
end

rawparams = odeTransformParams(params);

[tout, yout] = ode45(@(t, y) delayode(t, y, rawparams, p, t_obs), [t_obs], rawparams(5));
ll = gaussianLogLikelihood(yout, m, mvar + rawparams(6).^2);
if includeprior,
  ll = ll + odeLogPrior(rawparams);
end

if doplot,
  [tout, yout] = ode45(@(t, y) delayode(t, y, rawparams, p, t_obs), [min(t_obs), max(t_obs)], rawparams(5));
  %subplot(2, 1, 1);
  %plot(t_obs-min(t_obs), p, 'rx');
  %v = axis;
  %axis([0, 1300, 0, v(4)])
  %subplot(2, 1, 2);
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

grad = beta0 + beta*evalSpline(p, t_obs, t-Delta) - alpha * m;


function val = evalSpline(p, t_obs, t),

if t < t_obs(1),
  val = fnval(p.spline, log(p.timeshift + 0));
else
  val = fnval(p.spline, log(p.timeshift + t - t_obs(1)));
end


function ll = gaussianLogLikelihood(x, y, vars),

ll = sum(-.5*numel(x)*log(2*pi) - .5*sum(log(vars)) ...
	 -.5*sum((x-y).^2 ./ vars));


function ll = odeLogPrior(params),

bounds = [0, 1;
          1e-6, 1;
          1e-6, log(2);
          0, 299;
          0, 10;
          0, 10];

beta0 = params(1);
beta2 = params(2);
%alpha = log(2)/params(3);
alpha = params(3);
Delta = params(4);
m0 = params(5);
sigma_n = params(6);

ll = lnlogitnormalpdf(beta0, 0, 2, 0, 1) + ...
     lnlogitnormalpdf(beta2, 0, 2, 1e-6, 1) + ...
     lnlogitnormalpdf(alpha, 0, 2, 1e-6, log(2)) + ...
     lnlogitnormalpdf(Delta, -2, 2, 0, 300) + ...
     lnlogitnormalpdf(m0, 0, 2, 0, 10) + ...
     lnlogitnormalpdf(sigma_n, 0, 2, 0, 10);
