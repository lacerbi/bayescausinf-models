# Models of visuo-vestibular Bayesian causal inference

The file `bisensory_log_likelihood.m` computes the log-likelihood of a dataset `data` with model parameters `params`.

These are Bayesian causal inference models of visuo-vestibular perception applied to bisensory data for both localization (left/right) and unity judgement, stored in the file `bisensory_data.csv`. The CSV file contains 11 subjects (indexed by the first column).

Extract the desired dataset as follows from the CSV file:
```
id = 1; % Subject dataset (1-11)
data = csvread('bisensory_data.csv');   % Load data for all subjects
data_subj = data(data(:,1) == id,:);    % Get the target subject data
```

The considered model has 12 parameters:

1. Visual noise (low noise) - base standard deviation (log scale)
2. Visual noise (low noise) - Weber's fraction
3. Visual noise (medium noise) - base standard deviation (log scale)
4. Visual noise (medium noise) - Weber's fraction
5. Visual noise (high noise) - base standard deviation (log scale)
6. Visual noise (high noise) - Weber's fraction
7. Vestibular noise - base standard deviation (log scale)
8. Vestibular noise - Weber's fraction
9. Lapse rate
10. Gaussian prior mean
11. Gaussian prior standard deviation (log scale)
12. Probability of common cause (`p_common`)

Bounds for the parameters are provided below:
```
LB = [zeros(1,8), 0 -45 log(5) 0];
UB = [repmat([log(40),1], [1, 4]), 0.5 45 log(90) 1];
PLB = [repmat([log(2),0.02], [1, 4]), 0.02 -10 log(10) 0.1];
PUB = [repmat([log(20),0.3], [1, 4]), 0.2 10 log(45) 0.9];
```

The file `example_run.m` runs a maximum-likelihood optimization with BADS.

### Acknowledgments and references

The model file and data in this repository are obtained from [this repository](https://github.com/lacerbi/visvest-causinf), whose main reference is:

1. Acerbi\*, L., Dokka\*, K., Angelaki, D. E. & Ma, W. J. (2018). Bayesian Comparison of Explicit and Implicit Causal Inference Strategies in Multisensory Heading Perception, *PLoS Computational Biology* 14(7): e1006110. (\*equal contribution; [link](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006110))
