function ll = odeLikelihoodUT(p, m, mvar, t_obs, params, doplot),

if nargin < 6,
  doplot = 0;
end

[tout, yout] = ode45(@(t, y) delayode(t, y, params, p, t_obs), [t_obs], params(5));
ll = gaussianLogLikelihood(yout, m, mvar);

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
beta = params(2);
alpha = log(2)/params(3);
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
