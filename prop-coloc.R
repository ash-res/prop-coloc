#' Proportional colocalization test using genetic association summary data
#'
#' @description A proportional colocalization test that accounts for uncertainty in variant selection using one-sample summary data.
#'
#' @param b1 beta coefficients for trait 1
#' @param se1 standard errors for trait 1
#' @param b2 beta coefficients for trait 2
#' @param se2 standard errors for trait 2
#' @param n sample size
#' @param ld genetic variant correlation matrix
#' @param tau (optional) correlation between trait 1 and trait 2; if traits are measured in separate samples, then this should be set to 0 (default value is 0)
#' @param alpha (optional) nominal size of the conditional proportional colocalization test; if unspecified, p-value is reported to 2 decimal places
#' @param prune (optional) R^2 pruning threshold for variant correlations; if unspecified (default value is R^2 = 0.6)
#' @param J (optional) the top *J* variants for at least one trait are used to fit multivariable trait-variant linear models (default value is J=10)
#' @param figs (optional) logical: if \code{TRUE}, return plots of the univariable and fitted multivariable variant--trait associations of top *J* variants for each trait (default value is \code{FALSE})
#' @param traits (optional) a character vector of the two trait names
#' @param seed (optional) the seed for random sampling involved in computing conditional critical values for the prop-coloc-cond p-value (default value is 100)
#'
#' @details The method computes two different proportional colocalization tests. First, prop-coloc-cond is a conditional test based on lead variants and it accounts for uncertainty in variant selection. By contrast, prop-coloc-full tests the proportional colocalization hypothesis using a larger set of variant associations (by default, the top 10 strongest variant associations for each trait). Compared with prop-coloc-full, the prop-coloc-cond test is shown to achieve good type I error control, and be less sensitive to measurement errors in estimated genetic correlations. However, prop-coloc-cond and prop-coloc-full test slightly different null hypotheses; whereas prop-coloc-cond focuses on the evidence of lead variants, prop-coloc-full considers the evidence across more variants. We note that an optional input of "tau", which is the correlation between the two traits, is typically not provided in usual genetic association summary data. Therefore, sensitivity of analyses to the choice of tau is recommended.
#'
#' @return Output is a list containing:
#' 
#' \item{p_cond} p-value of prop-coloc-cond (to two decimal places), or (if \code{alpha} specified) logical \code{TRUE}/\code{FALSE} for whether the proportional colocalization null hypothesis is rejected (if \code{TRUE}) or not rejected (if \code{FALSE}) at the specified \code{alpha} significance threshold. Null hypothesis is that the beta-coefficients are proportional (this represents proportional colocalization), rejection indicates failure to proportionately colocalize.
#' \item{eta_cond} proportionality constant based on lead variants for each trait
#' \item{LM_cond} p-value of Lagrange Multiplier test based on lead variants for each trait. Null hypothesis is that the proportionality constant is zero, rejection indicates proportionality constant is non-zero. Colocalization is only meaningful if the proportionality constant is non-zero.
#' \item{p_full} p-value of prop-coloc-full test based on top *J* variants for each trait. Null hypothesis is that the beta-coefficients are proportional (this represents proportional colocalization), rejection indicates failure to colocalize.
#' \item{eta_full} proportionality constant based on top *J* variants for each trait
#' \item{LM_full} p-value of Lagrange Multiplier test based on top J variants for each trait. Null hypothesis is that the proportionality constant is zero, rejection indicates proportionality constant is non-zero. Colocalization is only meaningful if the proportionality constant is non-zero.
#' \item{Q} naive test statistic for proportional colocalization hypothesis from prop-coloc-cond
#' \item{alpha} nominal size of the conditional proportional colocalization test (if specified)
#' \item{convergence} whether uniroot function converged to find a conditional critical value for the prop-coloc-cond test
#' \item{top} two most relevant variants selected for the prop-coloc-cond test
#' \item{variants} set of variants used for the prop-coloc-full test
#' \item{fig_uni} plot of the univariable variant--trait associations of top *J* variants for each trait (if \code{figs} is specified and \code{TRUE})
#' \item{fig_multi} plot of the fitted multivariable variant--trait associations of top *J* variants for each trait (if \code{figs} is specified and \code{TRUE})
#'
#' @examples load("GLP1R.RData"); res <- prop.coloc(b1=thyroid$beta, se1=thyroid$se, b2=lung$beta, se2=lung$se, n=838, ld=ld)
#'
#' @references A frequentist test of proportional colocalization after selecting relevant genetic variants. Preprint.
#'
#' @export

prop.coloc <- function(b1,se1,b2,se2,n,ld,tau=NULL,alpha=NULL,prune=NULL,J=NULL,figs=NULL,traits=NULL,seed=100){
  if(missing(tau)){tau=0} else {tau=tau}
  if(missing(alpha)){alpha=seq(0.01,0.99,0.01)} else {alpha=alpha}
  if(missing(prune)){prune=0.6} else {prune=prune}
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
  
  # calculate t-statistics and re-order traits
  ts1 <- abs(b1)/se1
  ts2 <- abs(b2)/se2
  t2 <- (1*ifelse(sort(ts1,decreasing=TRUE)[1]>sort(ts2,decreasing=TRUE)[1],1,0))+(2*(1-ifelse(sort(ts1,decreasing=TRUE)[1]>sort(ts2,decreasing=TRUE)[1],1,0))) # trait 2 has the strongest lead variant
  t1 <- (1*ifelse(sort(ts1,decreasing=TRUE)[1]<sort(ts2,decreasing=TRUE)[1],1,0))+(2*(1-ifelse(sort(ts1,decreasing=TRUE)[1]<sort(ts2,decreasing=TRUE)[1],1,0)))
  if(t1==2){B1=b2; B2=b1; Se1=se2; Se2=se1; Ts1=ts2; Ts2=ts1}
  if(t1==2){b1=B1; b2=B2; se1=Se1; se2=Se2; ts1=Ts1; ts2=Ts2}
  if(t1==2){rm(B1,B2,Se1,Se2,Ts1,Ts2)}
  
  # pruning variants
  sel1 <- order(ts1,decreasing=TRUE)
  sel2 <- order(ts2,decreasing=TRUE)
  sel <- as.vector((t(cbind(sel1,sel2))))
  rm(sel1,sel2)
  colnames(ld) <- as.character(1:nrow(ld)); rownames(ld) <- as.character(1:nrow(ld)); snps <- as.character(1:nrow(ld))[sel]
  i=1
  while(i<=length(snps)){
    if(snps[i] %in% colnames(ld)){
      del <- which(as.vector((ld[which(colnames(ld)==snps[i]),-which(colnames(ld)==snps[i])]^2))>prune)
      del <- colnames(ld[-which(colnames(ld)==snps[i]),-which(colnames(ld)==snps[i])])[del] # snps to delete
      if(length(del)>0){
        snps <- snps[-which(snps %in% del)]
        b1 <- b1[-which(colnames(ld) %in% del)]; b2 <- b2[-which(colnames(ld) %in% del)]
        se1 <- se1[-which(colnames(ld) %in% del)]; se2 <- se2[-which(colnames(ld) %in% del)]
        ts1 <- ts1[-which(colnames(ld) %in% del)]; ts2 <- ts2[-which(colnames(ld) %in% del)]
        ld <- ld[-which(colnames(ld) %in% del),-which(colnames(ld) %in% del)]
      }
    }
    i=i+1
  }
  
  # consider J1 most relevant variants for each trait
  sel <- unique(c(order(ts1,decreasing=TRUE)[1:J1],order(ts2,decreasing=TRUE)[1:J1]))
  sel <- sort(sel)
  b1 <- b1[sel]; se1 <- se1[sel]; b2 <- b2[sel]; se2 <- se2[sel]; ld <- ld[sel,sel]; ts1 <- ts1[sel]; ts2 <- ts2[sel]; J <- length(b1)
  
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
  Sig11 <- solve(A1)*(1-as.numeric(t(B1)%*%solve(A1)%*%B1))
  Sig22 <- solve(A2)*(1-as.numeric(t(B2)%*%solve(A2)%*%B2))
  try(solve0 <- solve(sqrt.A1%*%t(sqrt.A2)))
  if(is.na(sum(solve0))){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample size n is too small.")}
  if(sum(c(diag(Sig11)<0,diag(Sig22)<0))>0){stop("Unable to compute the standard errors for the fitted multivariable regression model. \n This may be because the inputted sample size n is too small.")}
  Sig12 <- solve(sqrt.A1%*%t(sqrt.A2))*(tau-as.numeric(t(B1)%*%solve(sqrt.A1%*%t(sqrt.A2))%*%B2))
  Sig <- cbind(rbind(Sig11,t(Sig12)),rbind(Sig12,Sig22))
  gam1 <- as.vector(solve(A1)%*%B1)
  gam2 <- as.vector(solve(A2)%*%B2)
  gam <- c(gam1,gam2)
  
  # proportional colocalization test using all variants
  g <- function(eta){gam1-(gam2*eta)}
  Om <- function(eta){Sig11-(Sig12*eta)-(t(Sig12)*eta)+(Sig22*eta^2)}
  Q <- function(eta){as.numeric(t(g(eta))%*%solve(Om(eta))%*%g(eta))}
  init.val <- seq(-4,4,0.4)
  Q.init <- vector(,length=length(init.val))
  for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value}
  eta_full <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
  pval_full <- pchisq(n*Q(eta_full), df = J-1, lower.tail = F)
  
  # Lagrange multiplier test of eta=0 using all variants
  G = as.vector(-gam2)
  Om.0 <- Sig11; inv.Om.0 <- solve(Om.0); evec0 <- eigen(inv.Om.0)$vectors; eval0 <- eigen(inv.Om.0)$values
  sq.inv.Om.0 <- evec0%*%diag(sqrt(eval0))%*%t(evec0)
  LM_full <- n*as.numeric(t(t(G)%*%inv.Om.0%*%g(0))%*%solve(t(G)%*%inv.Om.0%*%G)%*%(t(G)%*%inv.Om.0%*%g(0)))
  LM_full <- pchisq(LM_full, df = 1, lower.tail = F)
  
  # Selecting the lead variants
  D <- diag(c(1/sqrt(diag(Sig11)),1/sqrt(diag(Sig22))))
  T0 <- as.vector(D%*%gam)*sqrt(n)
  T1 <- T0[1:J]; T2 <- T0[(J+1):(2*J)]
  j1 <- which.max(abs(T1))
  T2[j1] <- 0; j2 <- which.max(abs(T2))
  I. <- matrix(0,2,J); I.[1,j1] <- 1; I.[2,j2] <- 1
  
  # Proportional colocalization test using the lead variants only
  g. <- function(eta){as.vector(I.%*%(gam1-(gam2*eta)))}
  Om. <- function(eta){I.%*%(Sig11-(Sig12*eta)-(t(Sig12)*eta)+(Sig22*eta^2))%*%t(I.)}
  Q. <- function(eta){as.numeric(t(g.(eta))%*%solve(Om.(eta))%*%g.(eta))}
  Q.init <- vector(,length=length(init.val))
  for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q., method="Brent",lower=-1e2,upper=1e2)$value}
  eta.naive <- optim(init.val[which.min(Q.init)[[1]]], Q., method="Brent",lower=-1e2,upper=1e2)$par
  Q.naive <- Q.(eta.naive)
  
  # univariable and multivariable plots
  if(figs==TRUE){
    library(ggplot2)
    dat.beta <- data.frame(trait2 = b2, trait1 = b1)
    plots.beta <- ggplot(dat.beta, aes(x = trait2, y = trait1))
    p.beta <- plots.beta + geom_vline(xintercept = 0, color = "darkgrey", linetype = "twodash") + geom_hline(yintercept = 0, color = "darkgrey", linetype = "twodash") + geom_point(size=2) + geom_point(data=dat.beta[c(j1,j2),],pch=21, fill=NA, size=4, colour="forestgreen", stroke=1) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=11), plot.subtitle = element_text(hjust = 0.5), legend.key.size=unit(1, 'cm'), legend.text = element_text(size=12), legend.title = element_text(size=13)) + scale_x_continuous(limits=c(-max(abs(dat.beta)),max(abs(dat.beta)))) + scale_y_continuous(limits=c(-max(abs(dat.beta)),max(abs(dat.beta)))) + ylab(paste0(traits[t1])) + xlab(paste0(traits[t2])) + ggtitle("univariable (beta) associations")
    dat.gam <- data.frame(trait2 = gam2, trait1 = gam1)
    plots.gam <- ggplot(dat.gam, aes(x = trait2, y = trait1))
    p.gam <- plots.gam + geom_vline(xintercept = 0, color = "darkgrey", linetype = "twodash") + geom_hline(yintercept = 0, color = "darkgrey", linetype = "twodash") + geom_point(size=2) + geom_point(data=dat.gam[c(j1,j2),],pch=21, fill=NA, size=4, colour="forestgreen", stroke=1) + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=11), plot.subtitle = element_text(hjust = 0.5), legend.key.size=unit(1, 'cm'), legend.text = element_text(size=12), legend.title = element_text(size=13)) + scale_x_continuous(limits=c(-max(abs(dat.gam)),max(abs(dat.gam)))) + scale_y_continuous(limits=c(-max(abs(dat.gam)),max(abs(dat.gam)))) + ylab(paste0(traits[t1])) + xlab(paste0(traits[t2])) + ggtitle("multivariable (gamma) associations") + geom_abline(intercept = 0, slope = eta_full, color="firebrick1", linewidth=1) + geom_abline(intercept = 0, slope = eta.naive, color="forestgreen", linewidth=1)
  }

  # Conditional proportional colocalization test using the lead variants only
  G. <- as.vector(-I.%*%gam2)
  Om. <- Om.(eta.naive); inv.Om. <- solve(Om.); evec <- eigen(inv.Om.)$vectors; eval <- eigen(inv.Om.)$values
  sq.inv.Om. <- evec%*%diag(sqrt(eval))%*%t(evec)
  M. <- diag(2)-(sq.inv.Om.%*%G.%*%(1/(t(G.)%*%inv.Om.%*%G.))%*%t(G.)%*%sq.inv.Om.)
  C0 <- cbind(I.%*%(Sig11-(t(Sig12)*eta.naive)), I.%*%(Sig12-(Sig22*eta.naive)))
  C <- sq.inv.Om.%*%C0%*%t(D)
  
  # Lagrange Multipler test of eta=0
  Om.0 <- I.%*%Sig11%*%t(I.); inv.Om.0 <- solve(Om.0); evec0 <- eigen(inv.Om.0)$vectors; eval0 <- eigen(inv.Om.0)$values
  sq.inv.Om.0 <- evec0%*%diag(sqrt(eval0))%*%t(evec0)
  LM.0 <- n*as.numeric(t(t(G.)%*%inv.Om.0%*%g.(0))%*%solve(t(G.)%*%inv.Om.0%*%G.)%*%(t(G.)%*%inv.Om.0%*%g.(0)))
  LM.0 <- pchisq(LM.0, df = 1, lower.tail = F)
  
  # Conditional inference
  cond.test <- function(alpha){
  draws <- 1e4; K <- function(k){rnorm(2)}; K <- sapply(1:draws,K)
  ubar <- function(k){T0 - as.vector(t(C)%*%sq.inv.Om.%*%g.(eta.naive)*sqrt(n)) + as.vector(t(C)%*%K[,k])}
  ubar <- sapply(1:draws,ubar); ubar1 <- ubar[1:J,]; ubar2 <- ubar[(J+1):(2*J),]
  ubar2[j1,] <- 0
  A1 <- function(k){ifelse(which.max(abs(ubar1[,k]))==j1,1,0)*ifelse(which.max(abs(ubar2[,k]))==j2,1,0)}
  A1 <- sapply(1:draws,A1)
  A0 <- function(k){as.numeric(t(K[,k])%*%M.%*%K[,k])}
  A0 <- sapply(1:draws,A0)
  P <- function(w){(mean(ifelse(A0 <= w, 1,0)*A1)/mean(A1))-(1-alpha)}
  w0 <- NA
  try(w0 <- uniroot(P,c(0,500),extendInt="yes")$root, silent=T)
  if(is.na(w0)){w0 <- NA} else {w0 <- uniroot(P,c(0,500),extendInt="yes")$root}
  if(is.na(w0)){conv <- FALSE} else {conv <- TRUE}
  if(conv==FALSE){cond.res = FALSE}
  if(conv==TRUE){cond.res <- n*Q.naive>w0}
  return(c(cond.res,conv))
  }
  
  decimals <- function(number, places){format(round(number, places), nsmall = places)}
  
  if(length(alpha)==1){res <- cond.test(alpha); cond.res=res[1]; conv <- res[2]
    res.list <- list("p_cond"=cond.res, "eta_cond"=eta.naive, "LM_cond"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "LM_full"=LM_full, "Q"=Q.naive, "alpha"=alpha, "convergence"=conv, "top"=colnames(ld)[c(j1,j2)], "variants"=colnames(ld))
    Statistic <- c("Method", "eta", "p_LM", "p_coloc")
    res_cond <- c("prop-coloc-cond", decimals(res.list$eta_cond,2), decimals(res.list$LM_cond,2), res.list$p_cond)
    res_full <- c("prop-coloc-full", decimals(res.list$eta_full,2), decimals(res.list$LM_full,2), decimals(res.list$p_full,2))
    output.table <- data.frame(matrix(c(res_cond,res_full),nrow=2,byrow=TRUE))
    colnames(output.table) <- Statistic
    cat("\nProportional colocalization test\n")
    if (res.list$convergence == FALSE) { "\nConvergence failed for critical values for prop-coloc-cond test. Please interpret output with caution.\n" }
    cat("\n-----------------------------------\n")
    print(output.table, quote = F, row.names = FALSE, justify = "left")
    cat("-----------------------------------\n")
    cat("\np_coloc assesses proportional colocalization hypothesis (rejection indicates failure to colocalize).\n")
    if (res.list$p_cond==TRUE){cat("p_coloc = TRUE indicates proportional colocalization hypothesis is rejected at an alpha =",res.list$alpha,"significance level.\n")}
    if (res.list$p_cond==FALSE){cat("p_coloc = FALSE indicates proportional colocalization hypothesis is not rejected at an alpha =",res.list$alpha,"significance level.\n")}
    cat("p_LM assesses magnitude of the proportionality constant (rejection indicates proportionality constant is non-zero).\n")
    cat("Proportional colocalization is suggested when the colocalization hypothesis is not rejected and the proportionality constant is non-zero.\n")
  }
  
  if(length(alpha)>1){
    res <- sapply(alpha,cond.test); cond.res=res[1,]; conv <- res[2,]
      if(sum(cond.res)>0){
        res.list <- list("p_cond"=alpha[which(cond.res)[1]], "eta_cond"=eta.naive, "LM_cond"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "LM_full"=LM_full, "Q"=Q.naive, "convergence"=conv[which(cond.res)[1]], "top"=colnames(ld)[c(j1,j2)], "variants"=colnames(ld))
      }
      if(sum(cond.res)==0 & sum(conv)==length(alpha)){
        res.list <- list("p_cond"=max(alpha), "eta_cond"=eta.naive, "LM_cond"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "LM_full"=LM_full, "Q"=Q.naive, "convergence"=TRUE, "top"=colnames(ld)[c(j1,j2)], "variants"=colnames(ld))
      }
      if(sum(cond.res)==0 & sum(conv)!=length(alpha)){
        res.list <- list("p_cond"=FALSE, "eta_cond"=eta.naive, "LM_cond"=LM.0, "p_full"=pval_full, "eta_full"=eta_full, "LM_full"=LM_full, "Q"=Q.naive, "convergence"=FALSE, "top"=colnames(ld)[c(j1,j2)], "variants"=colnames(ld))
      }
    Statistic <- c("Method", "eta", "p_LM", "p_coloc")
    res_cond <- c("prop-coloc-cond", decimals(res.list$eta_cond,2), decimals(res.list$LM_cond,2), decimals(res.list$p_cond,2))
    res_full <- c("prop-coloc-full", decimals(res.list$eta_full,2), decimals(res.list$LM_full,2), decimals(res.list$p_full,2))
    output.table <- data.frame(matrix(c(res_cond,res_full),nrow=2,byrow=TRUE))
    colnames(output.table) <- Statistic
    cat("\nProportional colocalization test\n")
    if (res.list$convergence == FALSE) { "\nConvergence failed for critical values for prop-coloc-cond test. Please interpret output with caution.\n" }
    cat("\n-----------------------------------\n")
    print(output.table, quote = F, row.names = FALSE, justify = "left")
    cat("-----------------------------------\n")
    cat("\np_coloc assesses proportional colocalization hypothesis (rejection indicates failure to colocalize).\n")
    if (res.list$p_cond==TRUE){cat("p_coloc = TRUE indicates proportional colocalization hypothesis is rejected at an alpha =",res.list$alpha,"significance level.\n")}
    if (res.list$p_cond==FALSE){cat("p_coloc = FALSE indicates proportional colocalization hypothesis is not rejected at an alpha =",res.list$alpha,"significance level.\n")}
    cat("p_LM assesses magnitude of the proportionality constant (rejection indicates proportionality constant is non-zero).\n")
    cat("Proportional colocalization is suggested when the colocalization hypothesis is not rejected and the proportionality constant is non-zero.\n")
  }
  
  if(figs==TRUE){res.list$fig_uni=p.beta; res.list$fig_multi=p.gam}
  
  return(res.list)
}