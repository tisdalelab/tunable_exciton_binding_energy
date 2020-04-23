function [logprior] = makelogprior(lb,ub,nonlcon)
% This function will return logprior, a logical function which will be true
% if evaluated for a feasible walker.

% Inputs:
%   lb - lower bound vector for all parameters
%   ub - upper bound vector for all parameters
%   nonlcon - nonlinear constraints for the model of interest which takes a
%   walker as an input

% Output:
%   logprior - function handle which can be evaluated for a set of walkers.

% Check bounds are the same length
Nparams = length(lb);
if Nparams ~= length(ub)
    error('Upper and lower bounds are not the same length')
end

logprior = @(w) logical(true);

for i = 1:Nparams
    additional_constraint = @(w) w(i) < ub(i) & w(i) > lb(i);
    logprior = @(w) logprior(w) & additional_constraint(w);
end

% Add nonlinear constraints
if ~isempty(nonlcon)
    logprior = @(w) logprior(w) && nonlcon(w);
end

end

