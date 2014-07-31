rng(42);
N = 1000000;
for k=1:N,
  if rem(k, 1000) == 0,
    fprintf('%d/%d\n', k, N);
  end
  C = randn(202);
  C = C + C';
  [V, D] = eig(C);
end
