#"A frequentist test of proportional colocalization after selecting relevant genetic variants"#
#by Ashish Patel, John C. Whittaker, & Stephen Burgess#

* R code to run the method: prop-coloc.R
* R Markdown example: https://htmlpreview.github.io/?https://github.com/ash-res/prop-coloc/blob/main/prop-coloc.html

Interpreting key outputs:
 * full: p-value of prop-coloc-full test based on top J variants for each trait
 * * a low p-value suggests evidence *against* proportional colocalization when considering top J variants for each trait  
 * LM_full: p-value of LM test based on top J variants for each trait  
   a low p-value suggests evidence *for* a causal variant for trait 1 when considering top J variants for each trait  
 * cond: p-value of prop-coloc-cond (to two decimal places) if a specific significance level (alpha) is not inputted  
   a low p-value suggests evidence *against* proportional colocalization when considering lead variants for each trait
 * LM: p-value of LM test based on lead variants  
   a low p-value suggests evidence *for* a causal variant for trait 1 when considering lead variants for each trait

Inputs (*one* sample summary data):
  b1: beta coefficients for trait 1
  se1: standard errors for trait 1
  b2: beta coefficients for trait 2
  se2: standard errors for trait 2
  n: sample size
  ld: genetic variant correlation matrix
  tau (optional): correlation between trait 1 and trait 2; if unspecified, then set equal to 0
  alpha (optional): nominal size of the conditional proportional colocalisation test; otherwise reported to 2 decimal places
  prune (optional): R^2 pruning threshold for variant correlations; default: R^2 = 0.6
  J (optional): the top J variants for at least one trait are used to fit multivariable trait-variant linear models; default: J=10
  figs (optional): return plots of the fitted multivariable variant--trait associations of top J variants for each trait
  traits (optional): a character vector of the two trait names

Outputs:
  full: p-value of prop-coloc-full test based on top J variants for each trait  
  
  eta_full: eta (proportionality constant) estimate based on top J variants for each trait  
  
  LM_full: p-value of LM test based on top J variants for each trait  
  
  cond: p-value of prop-coloc-cond (to two decimal places) if alpha is not specified  
  ... or logical TRUE/FALSE for whether prop-coloc hypothesis is rejected (TRUE means the null is rejected)  
  
  eta: proportionality constant estimate based on lead variants  
  Q: GMM criterion based on lead variants
  LM: p-value of LM test based on lead variants
  convergence: whether uniroot function converged to find a conditional critical value for the prop-coloc-cond test
  top: two most relevant variants selected for the conditional proportional colocalization test
  variants: full set of variants used to fit multivariable linear variant--trait models
