# tesensitivity: A Stata package for assessing sensitivity to the unconfoundedness assumption

**Description**: A standard question in causal inference is to identify and estimate the effect of some treatment variable X on an outcome variable Y. 
A common assumption used to identify such effects is _unconfoundedness_, also known as selection on observables, conditional independence, ignorability, or exogenous selection. 
This assumption is not refutable, meaning that the data alone cannot tell us whether it is true. 
Nonetheless, researchers may wonder: How important is this assumption in their analyses? 
Put differently: How sensitive are results obtained under the conditional independence assumption? 

Masten and Poirier (2018) provide a set of theoretical results used to answer this question.
They define conditional _partial_ independence as any assumption weaker than full conditional independence. 
They specifically consider a class of assumptions called _conditional c-dependence_. 
This class measures relaxations of conditional independence by a single parameter c.
For any positive c, conditional independence only partially holds, and hence we cannot exactly learn the value of our treatment effect parameters, like ATE or ATT.
Instead, we only get bounds.
Masten and Poirier (2018) characterize these bounds as a function of c. 
Small values of c give narrow bounds while larger values of c give wider bounds.
Just how wide these bounds are, and hence how sensitive one's results are, depends on the data.

This repository contains a Stata module for estimating these bounds, and hence for checking the sensitivity of one's results to the conditional independence assumption. 
Masten and Poirier (2018) provide the formal definition, interpretation, and analysis of conditional c-dependence, and derive the theoretical identification results. 
The paper by Masten, Poirier, and Zhang (2020) explains the estimation strategy used by the Stata module.

**Authors**: This module was first written by Linqi Zhang (Boston College), and later extended by Paul Diegert (Duke), in collaboration with [Matt Masten](http://www.mattmasten.com) and [Alexandre Poirier](https://sites.google.com/site/alexpoirierecon/).

## Requirements

* Stata version 13

## Installation

Copy all files in `tesensitivityStataPackage` to your `Stata/ado/personal` directory. See the Stata help files for explanation of the syntax. 

## Subdirectories

* tesensitivityStataPackage - Contains all Stata package files
* vignette - Contains a vignette showing how to use the module. Specifically, it walks through the empirical illustration of LaLonde (1986) that is used in Masten, Poirier, and Zhang (2020).

## Troubleshooting

Please post problems or suggestions to the issue queue.

## References

LaLonde (1986) [Evaluating the Econometric Evaluations of Training Programs with Experimental Data](http://www.jstor.org/stable/1806062), _The American Economic Review_

Masten and Poirier (2018) [Identification of Treatment Effects under Conditional Partial Independence](https://mattmasten.github.io/assets/pdf/ECTA14481.pdff), _Econometrica_ ([Journal link](https://dx.doi.org/10.3982/ECTA14481))

Masten, Poirier, and Zhang (2020) [Assessing Sensitivity to Unconfoundedness: Estimation and Inference](https://arxiv.org/pdf/2012.15716), Working paper

## License

&copy; 2021 Linqi Zhang, Paul Diegert, Matt Masten, Alexandre Poirier

The contents of this repository are distributed under the MIT license. See file
`LICENSE` for details.
