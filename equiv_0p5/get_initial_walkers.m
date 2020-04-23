function [initial_walkers] = get_initial_walkers(lb,ub,logprior,Nwalkers)
% Initialize walkers within the supplied bounds

% Number of parameters
Nparams = length(lb);
Nwalkers_init = 2*Nwalkers;

if Nwalkers < 2*Nparams
    error('Use at least twice as many walkers as parameters.')
end

% Initialize
pinit = zeros(Nparams,Nwalkers_init);

% Randomly start walkers uniformly within the bounds
for n = 1:Nparams
        pinit(n,:) = unifrnd(lb(n),ub(n),[1,Nwalkers_init]);
end

% Eliminate walkers that do not satisfy the bounds. Will happen if there
% are nonlinear constraints included.
keep_vec = zeros(1,Nwalkers_init);
for n = 1:Nwalkers_init
    if logprior(pinit(:,n))
        keep_vec(n) = logical(true);
    end
end
pinit(:,~keep_vec) = [];

% Check that there are enough walker still
if size(pinit,2) < Nwalkers
    
    % If not enough then append the current walkers enough times to get to
    % the target number of walkers.
    while size(pinit,2) < Nwalkers
        pinit = [pinit,pinit];
        if size(pinit,2) > Nwalkers
            pinit(:,Nwalkers+1:end) = [];
            break
        end
    end
else
    
    % If there are more walkers than specified, only take as many as
    % desired.
    while size(pinit,2) > Nwalkers
        elim_ind = randi(size(pinit,2));
        pinit(:,elim_ind) = [];
    end
end

initial_walkers = pinit;

end

