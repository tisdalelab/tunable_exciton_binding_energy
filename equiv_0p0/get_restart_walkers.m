function [restart_walkers] = get_restart_walkers(pfull_prelim,logP_prelim)

% Number of walkers
Nwalkers = size(pfull_prelim,2);

% Find MLE
scale = 1.05;
[MLEval] = max(logP_prelim(2,:,end));
walker_inds_to_keep = logP_prelim(2,:,end) > scale*MLEval;

% Catch if there aren't enough good walkers
while sum(walker_inds_to_keep) < 3
    walker_inds_to_keep = logP_prelim(2,:,end) > scale*MLEval;
    scale = scale + 0.1;
end

% Take only those walkers within the close region to the MLE. Repeat that
% group of walkers to get back out to the correct number of walkers.
p_restart = pfull_prelim(:,walker_inds_to_keep,end);
while size(p_restart,2) < Nwalkers
    p_restart = [p_restart, p_restart];
    if size(p_restart,2) > Nwalkers
        p_restart(:,Nwalkers+1:end) = [];
        break
    end
end

restart_walkers = p_restart;

end

