% Written by Sam Winslow on March 9th, 2019
% Model function written and data collected by Katie Mauck

clear;
clc
home = pwd;

dataset = 5;
names = {'Eg','A','wL','wG','Ec','u','w2','E2','w3','E3','w4','E4','h','k','c1','c2','c3','B','bw','beep'};

%% User parameters

% Upper limit for eV to fit
maximum_eV_to_fit = 3.65;

% Bounds 
     %   Eg   A   wL    wG    Ec     u       w2    E2     w3    E3     w4    E4      h     k     c1    c2    c3    B      bw     beep    
ub =   [3.05  3   0.18  0.15  2.800  2.0     0.12  2.995  0.13  2.995  0.13  2.995   3.0   37    3.0   3.0   3.0   3.350  0.15   0.3  ];
lb =   [2.94  0.1 0.03  0.03  2.650  0.005   0.03  2.80   0.05  2.86   0.03  2.86    0.05  5     0.01  0.08  0.08  3.100  0.05   0.1  ];

% Error in absorbance data
abs_error = 0.02;

% MCMC specifications
Nwalkers_target = 200;
prelim_steps_per_walker = 5000;
production_steps_per_walker = 10000;
stepsize = 1.1;
burnin = 0.75;
confidence_level = 0.95;
parallel_spec = true;

% Output specs (true or false)
plot_prelim = false;
plot_production = false;

%% Get model and comprising functions
if dataset == 3
    [call_model,func_struct] = define_model_Q3();
elseif dataset == 4
    [call_model,func_struct] = define_model_Q4();
elseif dataset == 5
    [call_model,func_struct] = define_model_Q5();
end

%% load in all temperature points
A = csvread('Q5_AllT_A.csv');
T = csvread('Q5_AllT_T.csv');
x = csvread('Q5_AllT_x.csv');
X = 1240./x; % Convert wavelength to eV

N_T = length(T);

for i = 1:N_T
    
    %% Temperature specific operations
    Tspec = T(i);
    Tspec_str = num2str(Tspec);
    mkdir(Tspec_str)
    save_file_spec = ['Q5fitvals_linkedHx_',Tspec_str,'K'];

    %% Format data and select wavelengths to fit

    % Select x and y data for the lowest temperature for now (column 1)
    xdata = X(:,1);
    ydata = A(:,i);

    exclude_ind = find(xdata > maximum_eV_to_fit,1,'first');
    xdata = xdata(1:exclude_ind-1);
    ydata = ydata(1:exclude_ind-1);

    %% Setup loglikelihood and logprior functions

    % Make log-likelihood function as chi-squared statistics
    ystd = abs_error*ones(size(ydata));
    loglike = @(p) -sum( ( ydata - call_model(p,xdata) ).^2 ./ ystd.^2 );

    % Include bounds as the logprior
    nonlcon = @(p) p(8) < p(10) && p(10) < p(12) && p(12) < p(1);
    logprior = makelogprior(lb,ub,nonlcon);

    %% Setup preliminary MCMC run

    % Initialize walkers
    pinit = get_initial_walkers(lb,ub,logprior,Nwalkers_target);

    % Total number of steps to run for preliminary run
    totalprelimsteps = prelim_steps_per_walker*Nwalkers_target;

    % Run it
    tic
    [pfull_prelim,logP_prelim,Rej_prelim] = gwmcmc( pinit, { logprior loglike }, totalprelimsteps, 'burnin', 0, 'stepsize', stepsize ,'Parallel',logical(parallel_spec));
    prelim_time = toc;
    
    %% Plot preliminary run output

    if plot_prelim
        plot_MCMC_output(xdata,ydata,pfull_prelim,logP_prelim,call_model,func_struct,burnin,dataset,false)
    end

    %% Setup production run

    % Restart with best walkers
    p_restart = get_restart_walkers(pfull_prelim,logP_prelim);

    % Total number of steps to run for production run
    totalproductionsteps = production_steps_per_walker*Nwalkers_target;

    % Run it
    tic
    [pfull_final,logP_final,Rej_final] = gwmcmc( p_restart, { logprior loglike }, totalproductionsteps, 'burnin', 0, 'stepsize', stepsize ,'Parallel',logical(parallel_spec));
    production_time = toc;
    
    %% Plot production run output

    if plot_production
        plot_MCMC_output(xdata,ydata,pfull_final,logP_final,call_model,func_struct,burnin,dataset,true)
    end

    %% Extract statistics

    parametertable = parameter_statistics(pfull_final,logP_final,burnin,confidence_level,names);

    cd(Tspec_str)
    save(save_file_spec,'parametertable','pfull_prelim','logP_prelim','Rej_prelim','prelim_time','pfull_final','logP_final','Rej_final','production_time')
    cd(home)
end