function [parametertable] = parameter_statistics(pfull_final,logP,burnin,confidence_level,names)

% Select walkers past the burn-in
p = pfull_final(:,:,round(burnin*size(pfull_final,3)));

% Reshape to 2D array
p = p(:,:);

% Reorder parameters in each row lowest to highest
sortp = sort(p,2);

% Get confidence region
alpha = 1 - confidence_level;
leftCR = sortp(:,floor((alpha/2)*size(sortp,2)));
rightCR = sortp(:,ceil((1-alpha/2)*size(sortp,2)));

% Get mean parameter values
p_mean = mean(p,2);

% Get median parameter values
p_median = median(p,2);

% Get MLE parameters
[~,MLEind] = max(logP(2,:,end));
p_MLE = pfull_final(:,MLEind,end);

% Combine into a table
parametertable = array2table([leftCR,p_median,p_mean,p_MLE,rightCR],'VariableNames',{'LeftCR','Median','Mean','MLE','RightCR'},...
    'RowNames',names);
display(parametertable)
checkoutsiderange = (p_MLE > rightCR) | (p_MLE < leftCR);
if sum(checkoutsiderange) > 0
    warning('MLE is outside confidence region. Either not run long enough to be sampling a stationary distribution or a parameter is pushing up against a bound. Examine burn-in plot and corner plot and run again')
end

end

