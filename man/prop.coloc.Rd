% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prop-coloc.R
\name{prop.coloc}
\alias{prop.coloc}
\title{A frequentist test of proportional colocalization after selecting relevant genetic variants}
\usage{
prop.coloc(
  b1,
  se1,
  b2,
  se2,
  n,
  ld,
  tau = NULL,
  alpha = NULL,
  prune = NULL,
  J = NULL,
  figs = NULL,
  traits = NULL,
  seed = 100
)
}
\arguments{
\item{b1}{beta coefficients for trait 1}

\item{se1}{standard errors for trait 1}

\item{b2}{beta coefficients for trait 2}

\item{se2}{standard errors for trait 2}

\item{n}{sample size; a one-sample analysis is performed if only one sample size is entered. If two sample sizes are entered, then a two-sample analysis is performed under the assumption that genetic associations with traits 1 and 2 are measured from non-overlapping samples}

\item{ld}{genetic variant correlation matrix}

\item{tau}{(optional) correlation between trait 1 and trait 2; if traits are measured in separate samples, then this should be set to 0 (default value is 0)}

\item{alpha}{(optional) nominal size of the conditional proportional colocalization test; if unspecified, p-value is reported to 2 decimal places}

\item{prune}{(optional) R^2 pruning threshold for variant correlations; if unspecified (default value is R^2 = 0.6)}

\item{J}{(optional) the top J variants for at least one trait are used to fit multivariable trait-variant linear models (default value is J=10)}

\item{figs}{(optional) logical: if \code{TRUE}, return plots of the univariable and fitted multivariable variant--trait associations of top J variants for each trait (default value is \code{FALSE})}

\item{traits}{(optional) a character vector of the two trait names}

\item{seed}{(optional) the seed for random sampling involved in computing conditional critical values for the prop-coloc-cond p-value (default value is 100)}
}
\value{
Output is a list containing:
\itemize{
\item \code{p_cond} \cr
p-value of prop-coloc-cond test to two decimal places. Null hypothesis is that the beta-coefficients are proportional (this represents proportional colocalization), rejection indicates failure to proportionately colocalize.
\cr If \code{alpha} is specified, then \code{p_cond} is logical \code{TRUE} if the proportional colocalization null hypothesis is rejected or \code{FALSE} if the proportional colocalization null hypothesis is not rejected at the specified \code{alpha} significance threshold.
\item \code{eta_cond} \cr
proportionality constant based on lead variants for each trait
\item \code{LM_cond} \cr
p-value of Lagrange Multiplier test based on lead variants for each trait. Null hypothesis is that the proportionality constant is zero, rejection indicates proportionality constant is non-zero. Colocalization is only meaningful if the proportionality constant is non-zero.
\item \code{p_full} \cr
p-value of prop-coloc-full test based on top J variants for each trait. Null hypothesis is that the beta-coefficients are proportional (this represents proportional colocalization), rejection indicates failure to colocalize.
\item \code{eta_full} \cr
proportionality constant based on top J variants for each trait
\item \code{LM_full} \cr
p-value of Lagrange Multiplier test based on top J variants for each trait. Null hypothesis is that the proportionality constant is zero, rejection indicates proportionality constant is non-zero. Colocalization is only meaningful if the proportionality constant is non-zero.
\item \code{Q} \cr
naive test statistic for proportional colocalization hypothesis from prop-coloc-cond
\item \code{alpha} \cr
nominal size of the conditional proportional colocalization test (if specified)
\item \code{convergence} \cr
whether uniroot function converged to find a conditional critical value for the prop-coloc-cond test
\item \code{top} \cr
two most relevant variants selected for the prop-coloc-cond test
\item \code{variants} \cr
set of variants used for the prop-coloc-full test
\item \code{fig_uni} \cr
plot of the univariable variant--trait associations of top J variants for each trait (if \code{figs} is specified and \code{TRUE})
\item \code{fig_multi} \cr
plot of the fitted multivariable variant--trait associations of top J variants for each trait (if \code{figs} is specified and \code{TRUE})
}
}
\description{
A proportional colocalization test that accounts for uncertainty in variant selection using summary data.
}
\details{
The method computes two different proportional colocalization tests. Each test is a frequentist hypothesis test where the null hypothesis is colocalization (i.e. the genetic associations with traits 1 and 2 are proportional), and the alternative hypothesis is failure to colocalize (the genetic associations are not proportional). \cr
\cr
First, prop-coloc-cond is a conditional test based on lead variants for each trait that accounts for uncertainty in variant selection. By contrast, prop-coloc-full tests the proportional colocalization hypothesis using a larger set of variant associations (by default, the top 10 strongest variant associations for each trait).  \cr
\cr
Compared with prop-coloc-full, the prop-coloc-cond test is shown to achieve good type I error control, and be less sensitive to measurement errors in estimated genetic correlations. However, prop-coloc-cond and prop-coloc-full test slightly different null hypotheses; whereas prop-coloc-cond focuses on the evidence from lead variants, prop-coloc-full considers the evidence across more variants.  \cr
\cr
We note that an optional input of "tau", which is the correlation between the two traits, is typically not provided in usual genetic association summary data. Therefore, sensitivity of analyses to the choice of tau is recommended.  \cr
\cr
Both proportional colocalization tests require specification of the variant correlation matrix. This is a signed matrix, and the correlations must be provided with respect to the same effect alleles as the beta coefficients for each trait.
}
\examples{
res1 <- prop.coloc(b1=GLP1R$atrial_appendage$beta, se1=GLP1R$atrial_appendage$se, b2=GLP1R$pancreas$beta, se2=GLP1R$pancreas$se, ld=GLP1R$ld, n=838)
# here the LM test is rejected for both methods (p<0.01), so there is a non-zero slope in the scatterplot of genetic associations
# the proportionality test is rejected in the prop-coloc-full method (p<0.01), but not in the prop-coloc-cond method (p=0.10)
# hence there is evidence for failure to colocalize
# note that when alpha is not specified, the code takes much longer to run
res2 <- prop.coloc(b1=GLP1R$thyroid$beta, se1=GLP1R$thyroid$se, b2=GLP1R$lung$beta, se2=GLP1R$lung$se, ld=GLP1R$ld, n=838, alpha=0.05)
# here the LM test is not rejected, so there is no clear evidence for a non-zero slope in the scatterplot of genetic associations
# additionally, the proportionality test is rejected in both methods
# hence there is evidence for failure to colocalize
res3 <- prop.coloc(b1=GLP1R$atrial_appendage$beta, se1=GLP1R$atrial_appendage$se, b2=GLP1R$left_ventricle$beta, se2=GLP1R$left_ventricle$se, ld=GLP1R$ld, n=838, alpha=0.05, figs=TRUE)
# here the LM test is rejected for both methods (p<0.01), so there is a non-zero slope in the scatterplot of genetic associations
# the proportionality test is not rejected by either method: p>0.05 for prop-coloc-cond and p=0.72 for prop-coloc-full
# hence there is no evidence against colocalization
res3$fig_uni   # scatterplot of univariable associations (beta coefficients)
res3$fig_multi # scatterplot of multivariable associations

}
\references{
A frequentist test of proportional colocalization after selecting relevant genetic variants. Preprint. https://arxiv.org/abs/2402.12171
}
\author{
Stephen Burgess and Ashish Patel
}
