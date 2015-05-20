function r = merge_files(fnames),

% d = dir('ode_mcmc_2015-05-15_curves*.mat')
% r = merge_files({d.name})

res = cell(1, length(fnames));
for k=1:length(fnames),
  res{k} = load(fnames{k});
end

fields = fieldnames(res{1});

r = res{1};
for k=2:length(fnames),
  for l=1:length(fields),
    if iscell(res{k}.(fields{l})),
      I = ~cellfun(@isempty, res{k}.(fields{l}));
      r.(fields{l})(I) = res{k}.(fields{l})(I);
    else
      I = find(res{k}.(fields{l}));
      r.(fields{l})(I) = res{k}.(fields{l})(I);
    end
  end
end
