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

### Model parameters

The considered model has 12 parameters (elements of `params` array):

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
LB = [log(0.5)*ones(1,4), zeros(1,4), 0 -90 log(1) 0];
UB = [log(80)*ones(1,4), ones(1,4), 1 90 log(180) 1];
PLB = [log(1)*ones(1,4), 0.05*ones(1,4), 0.01 -5 log(4) 0.1];
PUB = [log(40)*ones(1,4), 0.5*ones(1,4), 0.2 5 log(90) 0.9];
```

### Example files

The file `example_run.m` runs a multi-start maximum-likelihood optimization with [BADS](https://github.com/acerbilab/bads).

### Acknowledgments and references

The model file and data in this repository are obtained from [this repository](https://github.com/lacerbi/visvest-causinf), whose main reference is:

1. Acerbi\*, L., Dokka\*, K., Angelaki, D. E. & Ma, W. J. (2018). Bayesian Comparison of Explicit and Implicit Causal Inference Strategies in Multisensory Heading Perception, *PLoS Computational Biology* 14(7): e1006110. (\*equal contribution; [link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006110))
