function [y, yvar] = unlog_bitseq_data(dlog),
  
% Logs are normally distributed
% ... recover mean in exp space.
y = exp(squeeze(dlog(:, 1, :) + dlog(:, 4, :)/2));

% Logs are normally distributed 
% ... recover variance in exp space.
yvar = squeeze((exp(dlog(:, 4, :))-1).*exp(2*dlog(:, 1, :) + dlog(:, 4, :)));
