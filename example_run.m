%RUN_OPTIMIZATION

id = 1; % Subject dataset (1-11)
data = csvread('bisensory_data.csv');   % Load data for all subjects
data_subj = data(data(:,1) == id,:);    % Get the target subject data

% Define parameter bounds
LB = [zeros(1,8), 0 -20 log(5) 0];
UB = [repmat([log(40),1], [1, 4]), 0.5 20 log(90) 1];
PLB = [repmat([log(2),0.01], [1, 4]), 0.02 -5 log(10) 0.1];
PUB = [repmat([log(20),0.5], [1, 4]), 0.2 5 log(45) 0.9];
x0 = 0.5*(PLB + PUB);

% Run maximum-likelihood estimation
options = bads('defaults');
options.UncertaintyHandling = false;
[x,nll] = bads(@(x) -bisensory_log_likelihood(x,data_one), x0, LB, UB, LB, UB, [], options);