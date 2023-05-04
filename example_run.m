%EXAMPLE_RUN Set up and run example maximum-likelihood estimation.

id = 1; % Subject dataset (1-11)
data = csvread('bisensory_data.csv');   % Load data for all subjects
data_subj = data(data(:,1) == id,:);    % Get the target subject data

% Define parameter bounds
LB = [0.5*ones(1,4), zeros(1,4), 0 -90 log(5) 0];
UB = [log(80)*ones(1,4), ones(1,4), 1 90 log(180) 1];
PLB = [log(1)*ones(1,4), 0.05*ones(1,4), 0.01 -5 log(4) 0.1];
PUB = [log(40)*ones(1,4), 0.5*ones(1,4), 0.2 5 log(90) 0.9];
Np = numel(PLB);    % Number of model parameters

%% Run maximum-likelihood estimation

% Target function: negative log-likelihood function for chosen dataset
fun = @(x) -bisensory_log_likelihood(x,data_subj);

% Optimizer options
options = bads('defaults');
options.UncertaintyHandling = false;
options.MaxFunEvals = 1;

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
    [x(i,:),nll(i)] = bads(fun, x0(i,:), LB, UB, LB, UB, [], options);
end