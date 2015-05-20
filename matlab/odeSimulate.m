function [yout, tout] = odeSimulate(p, t_p, t_obs, params, doplot),

if nargin < 5,
  doplot = 0;
end

%rawparams = odeTransformParams(params);
rawparams = params;

[tout, yout] = ode45(@(t, y) delayode(t, y, rawparams, p, t_p), [t_obs], rawparams{5});

if doplot,
  [tout, yout] = ode45(@(t, y) delayode(t, y, rawparams, p, t_p), [min(t_obs), max(t_obs)], rawparams{5});
  subplot(2, 1, 1);
  plot(t_p-min(t_p), p, 'rx');
  v = axis;
  axis([0, 1300, 0, v(4)])
  subplot(2, 1, 2);
  plot(tout-min(t_obs), yout);
  %hold on
  %plot(t_p-min(t_p), m, 'gx');
  %hold off
  v = axis;
  axis([0, 1300, 0, v(4)])
end


function grad = delayode(t, m, params, p, t_p),

beta0 = params{1};
beta2 = params{2};
beta = sqrt(beta2);
%alpha = log(2)/params(3);
if isscalar(params{3}),
  alpha = params{3};
else
  if t < params{3}(1),
    alpha = params{3}(2);
  else
    alpha = params{3}(3);
  end
end
Delta = params{4};

grad = beta0 + beta*joinTheDots(p, t_p, t-Delta) - alpha * m;


function val = joinTheDots(p, t_p, t),

if t < t_p(1),
  val = p(1);
else
  val = interp1(t_p, p, t, 'linear');
end
