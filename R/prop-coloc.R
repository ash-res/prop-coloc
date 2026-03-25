#' A frequentist test of proportional colocalization after selecting relevant genetic variants
#'
#' @description A proportional colocalization test that accounts for uncertainty in variant selection using summary data.
#'
#' @param b1 beta coefficients for trait 1
#' @param se1 standard errors for trait 1
#' @param b2 beta coefficients for trait 2
#' @param se2 standard errors for trait 2
#' @param n sample size; a one-sample analysis is performed if only one sample size is entered. If two sample sizes are entered, then a two-sample analysis is performed. Set tau=0 if traits 1 and 2 are measured from non-overlapping samples.
#' @param ld genetic variant correlation matrix
#' @param tau (optional) correlation between trait 1 and trait 2; if traits are measured in separate samples for a two-sample analysis, then this should be set to 0 (default value is 0)
#' @param pval (optional) nominal size of LM test for determining non-zero proportionality constant; if unspecified the default is pval = 0.05
#' @param alpha (optional) tested significance levels for determining p-value of the conditional proportional colocalization test; if unspecified the default is alpha = seq(0.01,0.99,0.01)
#' @param clump (optional) R^2 clumping threshold for variant correlations; if unspecified the default value is R^2 = 0.6
#' @param J (optional) the top J variants for at least one trait are used to fit multivariable trait-variant linear models (default value is J=10)
#' @param figs (optional) logical: if \code{TRUE}, return plots of the univariable and fitted multivariable variant--trait associations of top J variants for each trait (default value is \code{FALSE})
#' @param traits (optional) a character vector of the two trait names
#' @param seed (optional) the seed for random sampling involved in computing conditional critical values for the prop-coloc-cond p-value (default value is 100)
#'
#' @details The method computes two different proportional colocalization tests. Each test is a frequentist hypothesis test where the null hypothesis is colocalization (i.e. the genetic associations with traits 1 and 2 are proportional), and the alternative hypothesis is failure to colocalize (the genetic associations are not proportional). \cr
#' \cr
#' First, prop-coloc-cond is a conditional test based on lead variants for each trait that accounts for uncertainty in variant selection. By contrast, prop-coloc-full tests the proportional colocalization hypothesis using a larger set of variant associations (by default, the top 10 strongest variant associations for each trait).  \cr
#' \cr
#' Compared with prop-coloc-full, the prop-coloc-cond test is shown to achieve good type I error control, and be less sensitive to measurement errors in estimated genetic correlations. However, prop-coloc-cond and prop-coloc-full test slightly different null hypotheses; whereas prop-coloc-cond focuses on the evidence from lead variants, prop-coloc-full considers the evidence across more variants.  \cr
#' \cr
#' We note that an optional input of "tau", which is the correlation between the two traits, is typically not provided in usual genetic association summary data. Therefore, sensitivity of analyses to the choice of tau is recommended.  \cr
#' \cr
#' Both proportional colocalization tests require specification of the variant correlation matrix. This is a signed matrix, and the correlations must be provided with respect to the same effect alleles as the beta coefficients for each trait.
#'
#' @return Output is a list containing:
#' \itemize{
#'  \item \code{p_cond} \cr
#'  p-value of prop-coloc-cond test to two decimal places. Null hypothesis is that the beta-coefficients are proportional (this represents proportional colocalization), rejection indicates failure to proportionately colocalize.
#'  \cr If \code{alpha} is specified, then \code{p_cond} is logical \code{TRUE} if the proportional colocalization null hypothesis is rejected or \code{FALSE} if the proportional colocalization null hypothesis is not rejected at the specified \code{alpha} significance threshold.
#'  \item \code{eta_cond} \cr
#'  proportionality constant based on lead variants for each trait
#'  \item \code{p_eta} \cr
#'  p-value of Lagrange Multiplier test for a non-zero proportionality constant
#'  \item \code{p_full} \cr
#'  p-value of prop-coloc-full test based on top J variants for each trait. Null hypothesis is that the beta-coefficients are proportional (this represents proportional colocalization), rejection indicates failure to colocalize.
#'  \item \code{eta_full} \cr
#'  proportionality constant based on top J variants for each trait
#'  \item \code{p_naive} \cr
#'  naive p-value of prop-coloc-cond test statistic
#'  \item \code{alpha} \cr
#'  nominal size of the conditional proportional colocalization test (if specified)
#'  \item \code{convergence} \cr
#'  whether uniroot function converged to find a conditional critical value for the prop-coloc-cond test; a failure to converge may indicate a lack of support for a model of proportional colocalization
#'  \item \code{top} \cr
#'  two most relevant variants selected for the prop-coloc-cond test
#'  \item \code{variants} \cr
#'  set of variants used for the prop-coloc-full test
#'  \item \code{gamma} \cr
#'  the variant associations with each trait adjusted for all other variants included in a multivariable model
#'  \item \code{fig_uni} \cr
#'  plot of the univariable variant--trait associations of top J variants for each trait (if \code{figs} is specified and \code{TRUE})
#'  \item \code{fig_multi} \cr
#'  plot of the fitted multivariable variant--trait associations of top J variants for each trait (if \code{figs} is specified and \code{TRUE})
#' }
#'
#' @examples res1 <- prop.coloc(b1=GLP1R$stomach$beta, se1=GLP1R$stomach$se, b2=GLP1R$pancreas$beta, se2=GLP1R$pancreas$se, n=GLP1R$n$n_donors[c(8,7)], ld=GLP1R$ld, tau=GLP1R$trait_correlations[8,7],figs=TRUE,traits=c("stomach","pancreas"),clump=0.8)
#' # here the LM test is rejected (p<0.01; see \code{res1$p_eta}), so there is a non-zero slope in the scatterplot of genetic associations \code{res1$fig_multi}
#' # the proportionality test is rejected in the prop-coloc-full method (p<0.01; see \code{res1$p_full}), and in the prop-coloc-cond method (p=0.01)
#' # hence there is evidence for failure to colocalize
#' res2 <- prop.coloc(b1=GLP1R$stomach$beta, se1=GLP1R$stomach$se, b2=GLP1R$left_ventricle$beta, se2=GLP1R$left_ventricle$se, n=GLP1R$n$n_donors[c(8,4)], ld=GLP1R$ld, tau=GLP1R$trait_correlations[8,4],figs=TRUE,traits=c("stomach","left_ventricle"),clump=0.8)
#' # here the LM test is not rejected, so there is no clear evidence for a non-zero slope in the scatterplot of genetic associations \code{res2$fig_multi}
#' # hence the proportional colocalization hypothesis cannot be assessed since trait 1 may not have a causal variant
#' res3 <- prop.coloc(b1=GLP1R$atrial_appendage$beta, se1=GLP1R$atrial_appendage$se, b2=GLP1R$left_ventricle$beta, se2=GLP1R$left_ventricle$se, n=GLP1R$n$n_donors[c(2,4)], ld=GLP1R$ld, tau=GLP1R$trait_correlations[2,4],figs=TRUE,traits=c("atrial appendage","left_ventricle"),clump=0.8)
#' # here the LM test is rejected for both methods (p<0.01), so there is a non-zero slope in the scatterplot of genetic associations
#' # the proportionality test is not rejected by the prop-coloc-cond test (p=0.38)
#' # hence there is no evidence against proportional colocalization
#' \code{res3$fig_uni}   # scatterplot of univariable associations (beta coefficients)
#' \code{res3$fig_multi} # scatterplot of multivariable associations
#'
#' @references A frequentist test of proportional colocalization after selecting relevant genetic variants. Preprint. https://arxiv.org/abs/2402.12171
#'
#' @author Stephen Burgess and Ashish Patel
#'
#' @export

prop.coloc <- function(b1,se1,b2,se2,n,ld,tau=NULL,pval=NULL,alpha=NULL,clump=NULL,J=NULL,figs=NULL,traits=NULL,seed=100){
  if(missing(tau)){tau=0} else {tau=tau}
  if(missing(pval)){pval=0.05} else {pval=pval}
  if(missing(alpha)){alpha=seq(0.01,0.99,0.01)} else {alpha=alpha}
  if(missing(clump)){clump=0.6} else {clump=clump}
  if(missing(J)){J1=10} else {J1=J}
  if(missing(figs)){figs=FALSE}
  if(missing(traits)){traits = c("trait 1", "trait 2")}
  J = length(b1)
  if( exists(".Random.seed") ) {
    old <- .Random.seed
    on.exit( { .Random.seed <<- old } )
  }
  if (!is.na(seed)) { set.seed(seed) }

  # check the summary data are of the same dimension
  if(length(unique(c(length(b1),length(se1),length(b2),length(se2),nrow(ld),ncol(ld))))!=1){stop("The number of variants in the beta-coefficients, standard errors, and LD matrix are not the same.\n Please ensure only the same set of variants are included in each of these inputs.")}

  # clump variants
  clump.fun <- function(b1, se1, b2, se2, ld, clump){
    b1 <- as.vector(b1); se1 <- as.vector(se1)
    b2 <- as.vector(b2); se2 <- as.vector(se2)
    
    if(any(se1 == 0) | any(se2 == 0))
      stop("Standard errors contain zero")
    if(nrow(ld) != ncol(ld))
      stop("LD matrix must be square")
    
    p <- nrow(ld)
    
    if(length(b1) != p | length(se1) != p |
       length(b2) != p | length(se2) != p)
      stop("Effect vectors and LD matrix dimensions do not match")
    
    ts.1 <- abs(b1) / se1
    ts.2 <- abs(b2) / se2
    
    available <- rep(TRUE, p)
    keep      <- rep(FALSE, p)
    
    best1 <- which.max(ifelse(available, ts.1, -Inf))
    best2 <- which.max(ifelse(available, ts.2, -Inf))
    
    if(ts.1[best1] >= ts.2[best2]){
      current.trait <- 1L
      next.trait    <- 2L
    } else {
      current.trait <- 2L
      next.trait    <- 1L
    }
    
    ts.current <- if(current.trait == 1L) ts.1 else ts.2
    
    repeat {
      idx <- which.max(ifelse(available, ts.current, -Inf))
      if(!available[idx]) break
      keep[idx]      <- TRUE
      available[idx] <- FALSE
      r2        <- (ld[idx, ])^2
      to_remove <- which(r2 >= clump & available)
      available[to_remove] <- FALSE
      if(!any(available)) break
      ts.current <- if(next.trait == 1L) ts.1 else ts.2
      tmp           <- current.trait
      current.trait <- next.trait
      next.trait    <- tmp
    }
    
    list(
      keep = which(keep),
      b1  = b1[keep],
      se1 = se1[keep],
      b2  = b2[keep],
      se2 = se2[keep],
      ld  = ld[keep, keep]
    )
  }
  clump.res <- clump.fun(b1=b1, se1=se1, b2=b2, se2=se2, ld=ld, clump=clump)
  b1  <- clump.res$b1;  b2  <- clump.res$b2
  se1 <- clump.res$se1; se2 <- clump.res$se2
  ld  <- clump.res$ld
  ts1 <- abs(b1)/se1
  ts2 <- abs(b2)/se2
  rm(clump.res, clump.fun)

  # consider at most J1 relevant variants for each trait
  if(length(b1)<=J1){J1 <- length(b1)}
  sel <- unique(c(order(ts1,decreasing=TRUE)[1:J1],order(ts2,decreasing=TRUE)[1:J1]))
  sel <- sort(sel)
  b1 <- b1[sel]; se1 <- se1[sel]; b2 <- b2[sel]; se2 <- se2[sel]; ld <- ld[sel,sel]; ts1 <- ts1[sel]; ts2 <- ts2[sel]; J <- length(b1)

  ### ONE-SAMPLE ANALYSIS
  if(length(n)==1){
  # quantities of interest: trait 1
  a1 <- 1/((n*se1^2)+b1^2)
  A1 <- (sqrt(a1)%*%t(sqrt(a1)))*ld
  B1 <- a1*b1
  evec1 <- eigen(A1)$vectors; eval1 <- eigen(A1)$values
  sqrt.A1 <- evec1%*%diag(sqrt(eval1))%*%t(evec1)

  # quantities of interest: trait 2
  a2 <- 1/((n*se2^2)+b2^2)
  A2 <- (sqrt(a2)%*%t(sqrt(a2)))*ld
  B2 <- a2*b2
  evec2 <- eigen(A2)$vectors; eval2 <- eigen(A2)$values
  sqrt.A2 <- evec2%*%diag(sqrt(eval2))%*%t(evec2)
  
  # multivariable quantities
  try(solve0 <- solve(A1))
  if(is.na(sum(solve0))){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample sizes n are too small.")}
  try(solve0 <- solve(A2))
  if(is.na(sum(solve0))){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample sizes n are too small.")}
  try(solve0 <- solve(sqrt.A1%*%t(sqrt.A2)))
  if(is.na(sum(solve0))){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample size n is too small.")}
  Sig11 <- solve(A1)*(1-as.numeric(t(B1)%*%solve(A1)%*%B1))/n
  Sig22 <- solve(A2)*(1-as.numeric(t(B2)%*%solve(A2)%*%B2))/n
  if(sum(c(diag(Sig11)<0,diag(Sig22)<0))>0){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample size n is too small.")}
  if(tau == 0){Sig12 <- matrix(0,length(b1),length(b1))}
  if(tau != 0){Sig12 <- solve(sqrt.A1%*%t(sqrt.A2))*(tau-as.numeric(t(B1)%*%solve(sqrt.A1%*%t(sqrt.A2))%*%B2))/n}
  Sig <- cbind(rbind(Sig11,t(Sig12)),rbind(Sig12,Sig22))
  gam1 <- as.vector(solve(A1)%*%B1)
  gam2 <- as.vector(solve(A2)%*%B2)
  gam <- c(gam1,gam2)
  }
  
  if(length(n)==2){
  # quantities of interest: trait 1
  a1 <- 1/((n[1]*se1^2)+b1^2) # note: two-sample
  A1 <- (sqrt(a1)%*%t(sqrt(a1)))*ld
  B1 <- a1*b1
  evec1 <- eigen(A1)$vectors; eval1 <- eigen(A1)$values
  sqrt.A1 <- evec1%*%diag(sqrt(eval1))%*%t(evec1)
    
  # quantities of interest: trait 2
  a2 <- 1/((n[2]*se2^2)+b2^2) # note: two-sample
  A2 <- (sqrt(a2)%*%t(sqrt(a2)))*ld
  B2 <- a2*b2
  evec2 <- eigen(A2)$vectors; eval2 <- eigen(A2)$values
  sqrt.A2 <- evec2%*%diag(sqrt(eval2))%*%t(evec2)
  
  # multivariable quantities
  try(solve0 <- solve(A1))
  if(is.na(sum(solve0))){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample sizes n are too small.")}
  try(solve0 <- solve(A2))
  if(is.na(sum(solve0))){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample sizes n are too small.")}
  try(solve0 <- solve(sqrt.A1%*%t(sqrt.A2)))
  if(is.na(sum(solve0))){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample size n is too small.")}
  Sig11 <- solve(A1)*(1-as.numeric(t(B1)%*%solve(A1)%*%B1))/n[1]
  Sig22 <- solve(A2)*(1-as.numeric(t(B2)%*%solve(A2)%*%B2))/n[2]
  if(sum(c(diag(Sig11)<0,diag(Sig22)<0))>0){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample sizes n are too small.")}
  if(tau == 0){Sig12 <- matrix(0,length(b1),length(b1))}
  if(tau != 0){Sig12 <- solve(sqrt.A1%*%t(sqrt.A2))*(tau-as.numeric(t(B1)%*%solve(sqrt.A1%*%t(sqrt.A2))%*%B2))/sqrt(n[1]*n[2])}
  gam1 <- as.vector(solve(A1)%*%B1)
  gam2 <- as.vector(solve(A2)%*%B2)
  gam <- c(gam1,gam2)
  }

  # proportional colocalization test using all variants
  g <- function(eta){gam1-(gam2*eta)}
  Om <- function(eta){Sig11-(Sig12*eta)-(t(Sig12)*eta)+(Sig22*eta^2)}
  Q <- function(eta){as.numeric(t(g(eta))%*%solve(Om(eta))%*%g(eta))}
  init.val <- seq(-4,4,0.4)
  Q.init <- vector(,length=length(init.val))
  for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value}
  eta_full <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
  pval_full <- pchisq(Q(eta_full), df = J-1, lower.tail = F)

  # Selecting the lead variants
  G = as.vector(-gam2)
  D <- diag(c(1/sqrt(diag(Sig11)),1/sqrt(diag(Sig22))))
  T0 <- as.vector(D%*%gam)
  T1 <- T0[1:J]; T2 <- T0[(J+1):(2*J)]
  j2 <- which.max(abs(T2))
  T1[j2] <- 0; j1 <- which.max(abs(T1))
  I. <- matrix(0,2,J); I.[1,j1] <- 1; I.[2,j2] <- 1
  
  # Proportional colocalization test using the lead variants only
  g. <- function(eta){as.vector(I.%*%(gam1-(gam2*eta)))}
  Om. <- function(eta){I.%*%(Sig11-(Sig12*eta)-(t(Sig12)*eta)+(Sig22*eta^2))%*%t(I.)}
  Q. <- function(eta){as.numeric(t(g.(eta))%*%solve(Om.(eta))%*%g.(eta))}
  Q.init <- vector(,length=length(init.val))
  for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q., method="Brent",lower=-1e2,upper=1e2)$value}
  eta.naive <- optim(init.val[which.min(Q.init)[[1]]], Q., method="Brent",lower=-1e2,upper=1e2)$par
  Q.naive <- Q.(eta.naive)
  pval_naive <- pchisq(Q.naive, df = 2-1, lower.tail = F)
  
  # Lagrange Multipler test of eta=0
  G. <- as.vector(-I.%*%gam2)
  Om.0 <- I.%*%Sig11%*%t(I.); inv.Om.0 <- solve(Om.0)
  LM.0 <- as.numeric(t(t(G.)%*%inv.Om.0%*%g.(0))%*%solve(t(G.)%*%inv.Om.0%*%G.)%*%(t(G.)%*%inv.Om.0%*%g.(0)))
  LM.0 <- pchisq(LM.0, df = 1, lower.tail = F)
  
  # Check that lead variant for each trait is associated with that trait at p-value less than pval
  screen <- (LM.0 >= pval)
  
  # univariable and multivariable plots
  if(figs==TRUE){
    dat.beta <- data.frame(trait2 = b2, trait1 = b1)
    plots.beta <- ggplot2::ggplot(dat.beta, ggplot2::aes(x = trait2, y = trait1))
    p.beta <- plots.beta + ggplot2::geom_vline(xintercept = 0, color = "darkgrey", linetype = "twodash") + ggplot2::geom_hline(yintercept = 0, color = "darkgrey", linetype = "twodash") + ggplot2::geom_point(size=2) + ggplot2::geom_point(data=dat.beta[c(j1,j2),],pch=21, fill=NA, size=4, colour="forestgreen", stroke=1) + ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=11), plot.subtitle = ggplot2::element_text(hjust = 0.5), legend.key.size=ggplot2::unit(1, 'cm'), legend.text = ggplot2::element_text(size=12), legend.title = ggplot2::element_text(size=13)) + ggplot2::scale_x_continuous(limits=c(-max(abs(dat.beta)),max(abs(dat.beta)))) + ggplot2::scale_y_continuous(limits=c(-max(abs(dat.beta)),max(abs(dat.beta)))) + ggplot2::ylab(paste0(traits[1])) + ggplot2::xlab(paste0(traits[2])) + ggplot2::ggtitle("univariable (beta) associations")
    dat.gam <- data.frame(trait2 = gam2, trait1 = gam1)
    plots.gam <- ggplot2::ggplot(dat.gam, ggplot2::aes(x = trait2, y = trait1))
    p.gam <- plots.gam + ggplot2::geom_vline(xintercept = 0, color = "darkgrey", linetype = "twodash") + ggplot2::geom_hline(yintercept = 0, color = "darkgrey", linetype = "twodash") + ggplot2::geom_point(size=2) + ggplot2::geom_point(data=dat.gam[c(j1,j2),],pch=21, fill=NA, size=4, colour="forestgreen", stroke=1) + ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=11), plot.subtitle = ggplot2::element_text(hjust = 0.5), legend.key.size=ggplot2::unit(1, 'cm'), legend.text = ggplot2::element_text(size=12), legend.title = ggplot2::element_text(size=13)) + ggplot2::scale_x_continuous(limits=c(-max(abs(dat.gam)),max(abs(dat.gam)))) + ggplot2::scale_y_continuous(limits=c(-max(abs(dat.gam)),max(abs(dat.gam)))) + ggplot2::ylab(paste0(traits[1])) + ggplot2::xlab(paste0(traits[2])) + ggplot2::ggtitle("multivariable (gamma) associations") + ggplot2::geom_abline(intercept = 0, slope = eta_full, color="firebrick1", linewidth=1) + ggplot2::geom_abline(intercept = 0, slope = eta.naive, color="forestgreen", linewidth=1)
  }

  if(!screen){
  # Conditional proportional colocalization test using the lead variants only if they are associated with a trait at pval
  Om. <- Om.(eta.naive); inv.Om. <- solve(Om.); evec <- eigen(inv.Om.)$vectors; eval <- eigen(inv.Om.)$values
  sq.inv.Om. <- evec%*%diag(sqrt(eval))%*%t(evec)
  M. <- diag(2)-(sq.inv.Om.%*%G.%*%(1/(t(G.)%*%inv.Om.%*%G.))%*%t(G.)%*%sq.inv.Om.)
  C0 <- cbind(I.%*%(Sig11-(t(Sig12)*eta.naive)), I.%*%(Sig12-(Sig22*eta.naive)))
  C <- sq.inv.Om.%*%C0

  # Conditional inference
  # Pre-compute once
  draws <- 1e4
  K <- sapply(1:draws, function(k) rnorm(2))
  
  ubar <- sapply(1:draws, function(k) as.vector(gam) - as.vector(t(C)%*%sq.inv.Om.%*%g.(eta.naive)) + as.vector(t(C)%*%K[,k]))
  ubar1 <- ubar[1:J, ]
  ubar2 <- ubar[(J+1):(2*J), ]
  
  tbar <- sapply(1:draws, function(k) as.vector(D%*%ubar[,k]))
  tbar1 <- tbar[1:J, ]
  tbar2 <- tbar[(J+1):(2*J), ]
  tbar1[j2, ] <- 0
  
  lbar <- sapply(1:draws, function(k) as.numeric(t(t(I.%*%ubar2[,k])%*%inv.Om.0%*%(I.%*%ubar1[,k]))%*%solve(t(I.%*%ubar2[,k])%*%inv.Om.0%*%(I.%*%ubar2[,k]))%*%(t(I.%*%ubar2[,k])%*%inv.Om.0%*%(I.%*%ubar1[,k]))))
    A1 <- sapply(1:draws, function(k) {
      ifelse(which.max(abs(tbar1[,k])) == j1, 1, 0) *
      ifelse(which.max(abs(tbar2[,k])) == j2, 1, 0) * 
      ifelse(pchisq(lbar[k], df = 1, lower.tail = F) < pval, 1, 0)    
      })
    A0 <- sapply(1:draws, function(k) as.numeric(t(K[,k]) %*% M. %*% K[,k]))
  
  # Now alpha-specific part only
  cond.test <- function(alpha){
    P <- function(w) (mean(ifelse(A0 <= w, 1, 0) * A1) / mean(A1)) - (1 - alpha)
    w0 <- tryCatch(uniroot(P, c(0, 500), extendInt = "yes")$root, error = function(e) NA)
    conv <- !is.na(w0)
    cond.res <- if (!conv) FALSE else Q.naive > w0
    return(c(cond.res, conv))
  }
  }

  decimals <- function(number, places){format(round(number, places), nsmall = places)}

if(screen){
res.list <- list("p_cond"=NA, "eta_cond"=NA, "p_eta"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "p_naive"=NA, "alpha"=alpha, "convergence"=NA, "top"=NA, "variants"=colnames(ld), "gamma"=cbind(gam1,gam2))
Statistic <- c("Method", "eta", "p_coloc")
res_cond <- c("prop-coloc-cond", "NA", "NA")
res_full <- c("prop-coloc-full", decimals(res.list$eta_full,2), decimals(res.list$p_full,2))
output.table <- data.frame(matrix(c(res_cond,res_full),nrow=2,byrow=TRUE))
colnames(output.table) <- Statistic
cat("\nProportional colocalization test")
cat("\n-----------------------------------\n")
print(output.table, quote = F, row.names = FALSE, justify = "left")
cat("-----------------------------------\n")
cat("\np_coloc assesses proportional colocalization hypothesis (rejection indicates failure to colocalize).")
cat("\nProportional colocalization is suggested when the proportional colocalization hypothesis is not rejected.")
if (screen) { cat("\nLM test does not reject evidence for a non-zero proportionality constant at significance level ",pval)}
if (screen) { cat("\nOnly prop-coloc-full results are presented. Please interpret output with caution.")}
}
    
if(!screen){
  if(mean(A1)==0){
    res.list <- list("p_cond"=NA, "eta_cond"=eta.naive, "p_eta"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "p_naive"=pval_naive, "alpha"=alpha, "convergence"=NA, "top"=colnames(ld)[c(j1,j2)], "variants"=colnames(ld), "gamma"=cbind(gam1,gam2))
    Statistic <- c("Method", "eta", "p_coloc")
    res_cond <- c("prop-coloc-cond", "NA", "NA")
    res_full <- c("prop-coloc-full", decimals(res.list$eta_full,2), decimals(res.list$p_full,2))
    output.table <- data.frame(matrix(c(res_cond,res_full),nrow=2,byrow=TRUE))
    colnames(output.table) <- Statistic
    cat("\nProportional colocalization test")
    cat("\n-----------------------------------\n")
    print(output.table, quote = F, row.names = FALSE, justify = "left")
    cat("-----------------------------------\n")
    cat("\np_coloc assesses proportional colocalization hypothesis (rejection indicates failure to colocalize).")
    cat("\nProportional colocalization is suggested when the proportional colocalization hypothesis is not rejected.") 
    cat("\nA model of proportional colocalization based on lead variants is rejected. \nTherefore, prop-coloc-cond rejects the null hypothesis of proportional colocalization.") 
  }
  
  if((length(alpha)==1) && (mean(A1)>0)){
    res <- cond.test(alpha); cond.res=res[1]; conv <- res[2]
    res.list <- list("p_cond"=cond.res, "eta_cond"=eta.naive, "p_eta"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "p_naive"=pval_naive, "alpha"=alpha, "convergence"=conv, "top"=colnames(ld)[c(j1,j2)], "variants"=colnames(ld), "gamma"=cbind(gam1,gam2))
    Statistic <- c("Method", "eta", "p_coloc")
    res_cond <- c("prop-coloc-cond", decimals(res.list$eta_cond,2), res.list$p_cond)
    res_full <- c("prop-coloc-full", decimals(res.list$eta_full,2), decimals(res.list$p_full,2))
    output.table <- data.frame(matrix(c(res_cond,res_full),nrow=2,byrow=TRUE))
    colnames(output.table) <- Statistic
    cat("\nProportional colocalization test")
    cat("\n-----------------------------------\n")
    print(output.table, quote = F, row.names = FALSE, justify = "left")
    cat("-----------------------------------\n")
    cat("\np_coloc assesses proportional colocalization hypothesis (rejection indicates failure to colocalize).")
    cat("\nProportional colocalization is suggested when the proportional colocalization hypothesis is not rejected.")
    if (res.list$convergence == FALSE) { cat("\nConvergence failed for critical values for prop-coloc-cond test. Please interpret output with caution.")}
    if (res.list$p_cond==TRUE){cat("\np_coloc = TRUE indicates proportional colocalization hypothesis is rejected at an alpha =",res.list$alpha,"significance level.")}
    if (res.list$p_cond==FALSE){cat("\np_coloc = FALSE indicates proportional colocalization hypothesis is not rejected at an alpha =",res.list$alpha,"significance level.")}
  }

  if((length(alpha)>1) && (mean(A1)>0)){
    res <- sapply(alpha,cond.test); cond.res=res[1,]; conv <- res[2,]
      if(sum(cond.res)>0){
        res.list <- list("p_cond"=alpha[which(cond.res)[1]], "eta_cond"=eta.naive, "p_eta"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "p_naive"=pval_naive, "convergence"=conv[which(cond.res)[1]], "top"=colnames(ld)[c(j1,j2)], "variants"=colnames(ld), "gamma"=cbind(gam1,gam2))
      }
      if(sum(cond.res)==0 & sum(conv)==length(alpha)){
        res.list <- list("p_cond"=max(alpha), "eta_cond"=eta.naive, "p_eta"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "p_naive"=pval_naive, "convergence"=TRUE, "top"=colnames(ld)[c(j1,j2)], "variants"=colnames(ld), "gamma"=cbind(gam1,gam2))
      }
      if(sum(cond.res)==0 & sum(conv)!=length(alpha)){
        res.list <- list("p_cond"=FALSE, "eta_cond"=eta.naive, "p_eta"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "p_naive"=pval_naive, "convergence"=FALSE, "top"=colnames(ld)[c(j1,j2)], "variants"=colnames(ld), "gamma"=cbind(gam1,gam2))
      }
    Statistic <- c("Method", "eta", "p_coloc")
    res_cond <- c("prop-coloc-cond", decimals(res.list$eta_cond,2), if(is.numeric(res.list$p_cond)) decimals(res.list$p_cond,2) else "NA")
    res_full <- c("prop-coloc-full", decimals(res.list$eta_full,2), decimals(res.list$p_full,2))
    output.table <- data.frame(matrix(c(res_cond,res_full),nrow=2,byrow=TRUE))
    colnames(output.table) <- Statistic
    cat("\nProportional colocalization test")
    cat("\n-----------------------------------\n")
    print(output.table, quote = F, row.names = FALSE, justify = "left")
    cat("-----------------------------------\n")
    cat("\np_coloc assesses proportional colocalization hypothesis (rejection indicates failure to colocalize).")
    cat("\nProportional colocalization is suggested when the proportional colocalization hypothesis is not rejected.")
    if (res.list$convergence == FALSE) { cat("\nConvergence failed for critical values for prop-coloc-cond test. Please interpret output with caution.")}
    if (is.numeric(res.list$p_cond) & res.list$p_cond < max(alpha)){
      cat("\np_coloc =", decimals(res.list$p_cond, 2), 
          "indicates proportional colocalization hypothesis is rejected at this significance level.")
    } else if (is.numeric(res.list$p_cond) & res.list$p_cond == max(alpha)){
      cat("\np_coloc =", decimals(res.list$p_cond, 2), 
          "indicates proportional colocalization hypothesis is not rejected at any tested significance level.")
    } else {
      cat("\nConvergence failed; p_cond could not be determined.")
    }
  }
}

  if(figs==TRUE){res.list$fig_uni=p.beta; res.list$fig_multi=p.gam}

  return(res.list)
}
