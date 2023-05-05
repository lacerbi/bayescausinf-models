%EXAMPLE_RUN Set up and run example maximum-likelihood estimation.

% Use parametric Gaussian prior or semiparametric prior?
gaussian_prior = true;

id = 11; % Subject dataset (1-11)
data = csvread('bisensory_data.csv');   % Load data for all subjects
data_subj = data(data(:,1) == id,:);    % Get the target subject data

%% Define parameter bounds
LB = [log(0.5)*ones(1,4), zeros(1,4), 0 0, -90 log(1)];
UB = [log(80)*ones(1,4), ones(1,4), 1 1, 90 log(180)];
PLB = [log(1)*ones(1,4), 0.05*ones(1,4), 0.01 0.1, -5 log(4)];
PUB = [log(40)*ones(1,4), 0.5*ones(1,4), 0.2 0.9, 5 log(90)];

if ~gaussian_prior
    s_pivot = [0,1,2.5,5,10,15,25,35,45,60,90];
    num_pivots = numel(s_pivot) - 1;
    LB_prior = log(1e-3)*ones(1,num_pivots);
    PLB_prior = max(log(1e-2),log(-diff(-0.5*(s_pivot/90).^2)));
    PUB_prior = max(log(1),log(-diff(-0.5*(s_pivot/10).^2)));
    UB_prior = max(log(10),log(-diff(-0.5*(s_pivot/4).^2)));    
    LB = [LB(1:11),LB_prior];
    PLB = [PLB(1:11),PLB_prior];
    PUB = [PUB(1:11),PUB_prior];
    UB = [UB(1:11),UB_prior];
end

Np = numel(PLB);    % Number of model parameters

%% Run maximum-likelihood estimation

% Target function: negative log-likelihood function for chosen dataset
fun = @(x) -bisensory_log_likelihood(x,data_subj);

% Optimizer options
options = bads('defaults');
options.UncertaintyHandling = false;

% Run multi-start optimization
Nstarts = 3;

x0 = NaN(Nstarts,Np);
x = NaN(Nstarts,Np);
nll = NaN(1,Nstarts);

for i = 1:Nstarts
    if i == 1
        x0(i,:) = 0.5*(PLB + PUB); % Fixed first point
    else
        x0(i,:) = rand(1, numel(PLB)).*(PUB - PLB) + PLB; % Randomize
    end    
    [x(i,:),nll(i)] = bads(fun, x0(i,:), LB, UB, PLB, PUB, [], options);
end

% Temporary best results
% id = 11; x = [];

[idx,~] = min(nll);
x_best = x(idx,:);

%% Run Bayesian inference on all subjects

for id = 1:11
    data_subj = data(data(:,1) == id,:);    % Get the target subject data

    % Target function: log-joint for chosen dataset, with uniform prior
    log_prior = -sum(log(UB - LB));
    fun = @(x) bisensory_log_likelihood(x,data_subj) + log_prior;

    x0v = 0.5*(PLB + PUB);
    
    vbmc_options = vbmc('defaults');
    vbmc_options.MaxFunEvals = 700; % Default (100 + 50*nvars)
    [vp{id},elbo(id),elbo_sd(id),exitflag(id),output{id}] = ...
        vbmc(fun,x0v,LB,UB,PLB,PUB,vbmc_options);
end