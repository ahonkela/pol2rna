function p = lnlogitnormalpdf(x, mu, sigma2, a, b),

p = -0.5 * log(2*pi*sigma2) ...
    - 0.5 * ((logit((x-a)./(b-a)) - mu).^2 ./ sigma2) ...
    + log((b-a) ./ ((x-a) .* (b-x)));


function y = logit(x),

y = log(x ./ (1-x));
