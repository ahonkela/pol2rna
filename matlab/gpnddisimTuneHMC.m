function [epsilon, scales, m] = gpnddisimTuneHMC(m, options, scales),

oldparams = modelExtractParam(m);

if m.uniformPriors,
  if nargin < 3,
    bounds = gpnddisimExtractParamTransformSettings(m);
    bounds = cat(1, bounds{:})';
    basescales = diff(bounds) / max(diff(bounds));
    EPS_SCHEDULE = [0.0000001 0.000001 0.00001 0.0001 0.001 0.003 0.005 0.01 0.03 0.05 0.07 0.1 0.3 0.5 1 3 5 10];
    options.scales = basescales;
    options.epsilon = EPS_SCHEDULE(1);
    goodeps = options.epsilon;
    for myeps = 1:length(EPS_SCHEDULE),
      fprintf('Running with epsilon=%f\n', options.epsilon);
      [burninHMCsamples, Ehist] = gpnddisimSampleHMC(m, options.verbose, 100, options);
      %keyboard;
      if mean(all(diff(burninHMCsamples) == 0, 2)) > 0.2,
        break;
      else
        goodeps = options.epsilon;
        options.epsilon = EPS_SCHEDULE(find(EPS_SCHEDULE==goodeps)+1);
      end
      oldparams = burninHMCsamples(end,:);
      m = modelExpandParam(m, burninHMCsamples(end,:));
      disp(oldparams)
    end

    options.epsilon = 0.01;
    goodscales = 1e-6 * ones(size(oldparams));
    g = gpnddisimLogLikeGradients(m);
    I = find(g);
    for k = 1:length(I),
      fprintf('Tuning variable %d/%d\n', I(k), length(oldparams));
      accrate = 1;
      options.scales = 1e-6 * ones(size(oldparams));
      myscale = 0.01*basescales(I(k));
      while accrate > 0.8 && myscale < 25,
        myscale = myscale * 2;
        options.scales(I(k)) = myscale;
        fprintf('Variable %d/%d: scale %f\n', I(k), length(oldparams), myscale);
        [burninHMCsamples, Ehist] = gpnddisimSampleHMC(m, options.verbose, 20, options);
        %keyboard;
        %disp(g_mean)
        accrate = 1 - mean(all(diff(burninHMCsamples) == 0, 2));
        oldparams = burninHMCsamples(end,:);
        m = modelExpandParam(m, burninHMCsamples(end,:));
        disp(oldparams)
      end
      goodscales(I(k)) = myscale;
    end

    options.epsilon = goodeps;
    goodscales = 1e-6 * ones(size(oldparams));
    g = gpnddisimLogLikeGradients(m);
    I = find(g);
    for k = 1:length(I),
      fprintf('Tuning variable %d/%d\n', I(k), length(oldparams));
      accrate = 1;
      options.scales = 1e-6 * ones(size(oldparams));
      myscale = 0.1*basescales(I(k));
      while accrate > 0.8 && myscale < 25,
        myscale = myscale * 2;
        options.scales(I(k)) = myscale;
        fprintf('Variable %d/%d: scale %f\n', I(k), length(oldparams), myscale);
        [burninHMCsamples, Ehist] = gpnddisimSampleHMC(m, options.verbose, 20, options);
        %keyboard;
        %disp(g_mean)
        accrate = 1 - mean(all(diff(burninHMCsamples) == 0, 2));
        oldparams = burninHMCsamples(end,:);
        m = modelExpandParam(m, burninHMCsamples(end,:));
        disp(oldparams)
      end
      goodscales(I(k)) = myscale;
    end
  else
    goodscales = scales;
  end

  fprintf('Scales:\n');
  disp(goodscales)
else
  goodscales = ones(size(oldparams));
end

EPS_SCHEDULE = [0.00001 0.0001 0.001 0.003 0.005 0.01 0.03 0.05 0.07 0.1 0.3 0.5 1 3 5 10];
options.epsilon = EPS_SCHEDULE(1);
goodeps = options.epsilon;
options.scales = goodscales;
for myeps = 1:length(EPS_SCHEDULE),
  fprintf('Running with epsilon=%f\n', options.epsilon);
  [burninHMCsamples, Ehist] = gpnddisimSampleHMC(m, options.verbose, 100, options);
  %keyboard;
  rejrate = mean(all(diff(burninHMCsamples) == 0, 2));
  if rejrate > 0.2,
    fprintf('Too high rejection rate (%.1f%%), backtracking eps schedule...\n', 100*rejrate);
    options.epsilon = goodeps;
  else
    fprintf('Good rejection rate (%.1f%%)\n', 100*rejrate);
    goodeps = options.epsilon;
    options.epsilon = EPS_SCHEDULE(find(EPS_SCHEDULE==goodeps)+1);
  end
  oldparams = burninHMCsamples(end,:);
  m = modelExpandParam(m, burninHMCsamples(end,:));
  disp(oldparams)
end
epsilon = goodeps;
scales = goodscales;
