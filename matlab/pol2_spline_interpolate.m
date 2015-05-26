% POL2_SPLINE_INTERPOLATE
%
% Testing Matlab cubic spline smoothing with Pol-II data and
% using leave-one-out cross validation to tune the regularisation

% (c) 2015 Antti Honkela

d = load('ode_data_2015-05-07.mat');

N = size(d.dataVals1, 2);
K = 8;
p_test = linspace(0.8, 1, 21);
tshift_test = 5;

myfits = cell(K, N);
pol2fits = cell(N, 1);
errors = zeros(K, N);
meanerrors = zeros(length(p_test), length(tshift_test));
geneerrors = zeros(length(p_test), N);
for j=1:length(tshift_test),
  timevector = log(d.timevector - 295);
  %timevector = d.timevector - 295;
  for l=1:length(p_test),
    fprintf('j=%d, doing test %d/%d\n', j, l, length(p_test));
    for k=1:K,
      I = setdiff(1:10, k+1);
  
      for n=1:N,
        myfit = csaps(timevector(I), d.dataVals1(I,n), p_test(l));
        errors(k, n) = fnval(myfit, timevector(k+1)) - d.dataVals1(k+1,n);
      end
    end
    meanerrors(l,j) = mean(mean(errors'.^2));
    geneerrors(l,:) = mean(errors.^2);
  end
end
[~, J] = min(meanerrors);
for n=1:N,
  pol2fits{n} = csaps(timevector, d.dataVals1(:,n), p_test(J));
end

splinevars = geneerrors(J,:) ./ var(d.dataVals1);
save('ode_spline_vars.mat', 'splinevars');
