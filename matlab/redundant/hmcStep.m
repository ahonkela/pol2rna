function [x, E] = hmcStep(findE, x, options, gradE, varargin),

L = options.maxit;

Tau = options.tau;
epsilon = options.epsilon;

g = feval(gradE, x, varargin{:}); % set gradient using initial x 
E = feval(findE, x, varargin{:}); % set objective function too 

p = randn(size(x));
H = p * p' / 2 + E ; % evaluate H(x,p) 
xnew = x ; gnew = g ; 
try,
  for tau = 1:Tau % make Tau 'leapfrog' steps 
    p = p - epsilon * gnew / 2 ; % make half-step in p 
    xnew = xnew + epsilon * p ; % make step in x 
    gnew = feval(gradE, xnew, varargin{:}); % find new gradient 
    p = p - epsilon * gnew / 2 ; % make half-step in p 
  end 
  Enew = feval(findE, xnew, varargin{:}) ; % find new value of H 
  Hnew = p * p' / 2 + Enew ; 
  dH = Hnew - H ; % Decide whether to accept 
  if options.verbose,
    fprintf('Step, threshold: %.4f\n', exp(-dH));
  end
  if ( dH < 0 ) accept = 1 ; 
  elseif ( rand() < exp(-dH) ) accept = 1 ; 
  else accept = 0 ; 
  end 
  if ( accept ) 
    if options.verbose,
      fprintf('accepted.\n');
      disp(x);
    end
    g = gnew ; x = xnew ; E = Enew ; 
  end 
catch me,
  rethrow(me);
end
