function [mean, vars] = gpnddisimSamplePredictions(model, paramsamples, t_pred, N_samples),

% t_pred = (((0:100)/100*sqrt(1280)).^2 + 300)';

if nargin < 4,
  N_samples = 10;
end

r = zeros(size(paramsamples, 1)*N_samples, 2*length(t_pred));
baseind = 0;
for k=1:size(paramsamples, 1),
  m = modelExpandParam(model, paramsamples(k, :));
  [~, mu, C] = gpnddisimPredict(m, t_pred, 1, 0);
  C = 0.5*(C+C');
  eigdone = 0;
  jitter = 0;
  jitterincr = 1e-14;
  % Work around a Matlab EIG bug that sometimes causes it to fail with
  % message "EIG did not converge"
  while ~eigdone,
    try,
      [V,D]=eig(C + jitter * eye(size(C)));
      eigdone = 1;
    catch
      fprintf('EIG did not converge with jitter %g, trying again...\n', jitter);
      %tmp_name = tempname('/home/fs/ahonkela/mlprojects/pol2rnaseq/matlab');
      %save(tmp_name, 'C');
      jitter = jitter + jitterincr;
      jitterincr = 2*jitterincr;
    end
  end
  D = diag(diag(D) - jitter);
  if min(diag(D)) < -1e-8,
    fprintf('min(D) = %g\n', min(diag(D)));
  end
  if jitter > 0,
    r(baseind+1:baseind+N_samples, :) = mvnrnd(mu', C + jitter * eye(size(C)), N_samples);
  else
    eigdone = 0;
    jitter = 0;
    jitterincr = 1e-14;
    while ~eigdone,
      try,
        r(baseind+1:baseind+N_samples, :) = mvnrnd(mu', C + ((jitter-1.1*min(0, min(diag(D)))) * eye(size(C))), N_samples);
        eigdone = 1;
      catch
        fprintf('EIG did not converge with jitter %g, trying again...\n', jitter);
        jitter = jitter + jitterincr;
        jitterincr = 2*jitterincr;
      end
    end
  end
  baseind = baseind + N_samples;
end
mean = r;
vars = [];
