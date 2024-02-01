## A frequentist test of proportional colocalization after selecting relevant genetic variants ##
#### by Ashish Patel, John C. Whittaker, & Stephen Burgess #### 


##### Installation
* devtools::install_github("ash-res/prop-coloc")

##### R Markdown example
* https://htmlpreview.github.io/?https://github.com/ash-res/prop-coloc/blob/main/prop-coloc.html

##### Interpreting key outputs
 * p_cond: p-value of prop-coloc-cond (to two decimal places) if a specific significance level (alpha) is not specified
     * a low p-value suggests evidence *against* proportional colocalization when considering the lead variants for each trait
 * LM_cond: p-value of LM test based on lead variants
     * a low p-value suggests evidence *for* a causal variant for trait 1 when considering lead variants for each trait
 * p_full: p-value of prop-coloc-full test based on top J variants for each trait
     * a low p-value suggests evidence *against* proportional colocalization when considering top J variants for each trait  
 * LM_full: p-value of LM test based on top J variants for each trait
     * a low p-value suggests evidence *for* a causal variant for trait 1 when considering top J variants for each trait

##### Inputs (genetic association summary data)
 * b1: beta coefficients for trait 1
 * se1 standard errors for trait 1
 * b2: beta coefficients for trait 2
 * se2: standard errors for trait 2
 * n: sample size
 * ld: genetic variant correlation matrix
 * tau (optional): correlation between trait 1 and trait 2; if traits are measured in separate samples, then this should be set to 0 (default value is 0)
 * alpha (optional): nominal size of the conditional proportional colocalization test; if unspecified, p-value is reported to 2 decimal places
 * prune (optional): R^2 pruning threshold for variant correlations; if unspecified (default value is R^2 = 0.6)
 * J (optional): the top J variants for at least one trait are used to fit multivariable trait-variant linear models (default value is J=10)
 * figs (optional): logical: if *TRUE*, return plots of the univariable and fitted multivariable variant--trait associations of top J variants for each trait (default value is *FALSE*)
 * traits (optional): a character vector of the two trait names
 * seed (optional): the seed for random sampling involved in computing conditional critical values for the prop-coloc-cond p-value (default value is 100)

##### Output is a list containing...
 * p_cond: p-value of prop-coloc-cond (to two decimal places), or (if *alpha* specified) logical *TRUE*/*FALSE* for whether the proportional colocalization null hypothesis is rejected (if *TRUE*) or not rejected (if *FALSE*) at the specified *alpha* significance threshold. Null hypothesis is that the beta-coefficients are proportional (this represents proportional colocalization), rejection indicates failure to proportionately colocalize.
 * eta_cond: proportionality constant based on lead variants for each trait
 * LM_cond: p-value of Lagrange Multiplier test based on lead variants for each trait. Null hypothesis is that the proportionality constant is zero, rejection indicates proportionality constant is non-zero. Colocalization is only meaningful if the proportionality constant is non-zero.
 * p_full: p-value of prop-coloc-full test based on top J variants for each trait. Null hypothesis is that the beta-coefficients are proportional (this represents proportional colocalization), rejection indicates failure to colocalize.
 * eta_full: proportionality constant based on top J variants for each trait
 * LM_full: p-value of Lagrange Multiplier test based on top J variants for each trait. Null hypothesis is that the proportionality constant is zero, rejection indicates proportionality constant is non-zero. Colocalization is only meaningful if the proportionality constant is non-zero.
 * Q: naive test statistic for proportional colocalization hypothesis from prop-coloc-cond
 * alpha: nominal size of the conditional proportional colocalization test (if specified)
 * convergence: whether uniroot function converged to find a conditional critical value for the prop-coloc-cond test
 * top: two most relevant variants selected for the prop-coloc-cond test
 * variants: set of variants used for the prop-coloc-full test
 * fig_uni: plot of the univariable variant--trait associations of top J variants for each trait (if *figs* is specified and *TRUE*)
 * fig_multi: plot of the fitted multivariable variant--trait associations of top J variants for each trait (if *figs* is specified and *TRUE*)

##### Example
 * res <- prop.coloc(b1=GLP1R$thyroid$beta, se1=GLP1R$thyroid$se, b2=GLP1R$lung$beta, se2=GLP1R$lung$se, n=838, ld=ld, alpha=0.05)
