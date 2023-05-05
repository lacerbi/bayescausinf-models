# Models of visuo-vestibular Bayesian causal inference

We provide here example files to easily run and test algorithms for fitting models of Bayesian causal inference in multisensory perception, here in the visuo-vestibular context (Acerbi et al., 2018).

## Description of the files

The file `bisensory_log_likelihood.m` computes the log-likelihood of a dataset `data` with model parameters `params`.

### Data

The data are bisensory (visuo-vestibular) datasets for both localization (left/right) and unity judgement, stored in the file `bisensory_data.csv`. The CSV file contains 11 subjects (indexed by the first column).

For a given row, the columns contain data for:
1. Subject id
2. Number of trial
3. Task type (1 unused, 2 Vestibular localisation, 3 Unity judgment)
4. Vestibular stimulus position (deg)
5. Visual stimulus position (deg)
6. Response (1 or -1 for right/left; 1 or 2 for yes/no unity)
7. Visual noise level (1 low, 2 med, 3 high)

Extract the desired dataset as follows from the CSV file:
```
id = 1; % Subject dataset (1-11)
data = csvread('bisensory_data.csv');   % Load data for all subjects
data_subj = data(data(:,1) == id,:);    % Get the target subject data
```
### Models

The model implemented is a Bayesian causal inference observer model of multisensory perception with stimulus-dependent noise, which induces non-Gaussian likelihoods for the observer's sensory measurements. The Bayesian observer's posterior is computed via numerical integration.

At the moment, `bisensory_log_likelihood.m` implements two different models for the observer's prior:
- The first model is a standard Gaussian prior (12 parameters in total).
- The second model is a semiparametric prior, defined on a grid of pivots (21 parameters in total).

The function distinguishes between the two models based on the number of elements of the `params` array.

#### Gaussian prior model (12 parameters)

The considered observer model has 12 parameters in total (elements of `params` array):

1. Noise base standard deviation, visual modality (low noise condition) - log scale
2. Noise base standard deviation, visual modality (medium noise condition) - log scale
3. Noise base standard deviation, visual modality (high noise condition) - log scale
4. Noise base standard deviation, vestibular modality - log scale
5. Noise Weber's fraction, visual modality (low noise condition)
6. Noise Weber's fraction, visual modality (medium noise condition)
7. Noise Weber's fraction, visual modality (high noise condition)
8. Noise Weber's fraction, vestibular modality
9. Lapse rate (probability of random response)
10. Probability of common cause (`p_common`)
11. Gaussian prior mean
12. Gaussian prior standard deviation (log scale)

Bounds for the parameters are provided below (these are the same as Acerbi et al., 2018):
```
LB = [log(0.5)*ones(1,4), zeros(1,4), 0 0, -90 log(1)];
UB = [log(80)*ones(1,4), ones(1,4), 1 1, 90 log(180)];
PLB = [log(1)*ones(1,4), 0.05*ones(1,4), 0.01 0.1, -5 log(4)];
PUB = [log(40)*ones(1,4), 0.5*ones(1,4), 0.2 0.9, 5 log(90)];
```

#### Semiparametric prior model (21 parameters)

This model has 21 parameters in total. The first 11 parameters have the same meaning as the basic model defined above (where the 11th parameter is the prior mean).
Then the model has 10 further parameters which encode differences in the log-density of the prior, evaluated on a pre-specified pivot grid.
The prior is assumed to be symmetric around the mean and monotonically decreasing from the mean (although it can be near-flat).

Code to define the parameters bounds is provided below, to add to the code previously provided for the Gaussian prior model:
```
% LB, UB, PLB, PUB have been defined before
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
```

### Example files

The file `example_run.m` runs a multi-start maximum-likelihood optimization with [BADS](https://github.com/acerbilab/bads).

### Acknowledgments and references

The model file and data in this repository are obtained from [this repository](https://github.com/lacerbi/visvest-causinf), whose main reference is:

1. Acerbi\*, L., Dokka\*, K., Angelaki, D. E. & Ma, W. J. (2018). Bayesian Comparison of Explicit and Implicit Causal Inference Strategies in Multisensory Heading Perception, *PLoS Computational Biology* 14(7): e1006110. (\*equal contribution; [link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006110))
