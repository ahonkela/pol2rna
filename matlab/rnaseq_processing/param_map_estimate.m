function theta_hat=param_map_estimate(negloglike,theta_start)

% usage: theta_hat = param_map_estimate(@(th) negloglike(th,y), theta)


theta_hat = fminsearch(@(th) negloglike(th),theta_start);
end
