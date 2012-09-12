function multipriorDisplay(model, spacing_count)

% MULTIPRIORDISPLAY Display a multiprior model.

if nargin > 1
  spacing = repmat(32, 1, spacing_count);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Multiprior model:\n')
modelDisplay(model.basemodel, spacing_count+2);
