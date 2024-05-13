############################################################
# Section 3: selective CIs in multiple testing
# helper scripts
############################################################
library(MASS)

gen_data <- function(tau) {
  ######################################
  # Generate Gaussian noise w specific
  # correlation structure using
  # Toeplitz matrix
  ######################################
  
  # Create data and fission split
  Sigma <- diag(rep(1, n))
  Y <- mvrnorm(1, mu, Sigma)
  Z <- mvrnorm(1, rep(0, n), Sigma)
  f_Y <- Y + tau*Z
  g_Y <- Y - (1/tau)*Z
  
  # Also return p-values for fissioned and non-fissioned data
  pvals <- 1 - pnorm(Y)
  pvals_fission <- 1 - pnorm(f_Y, sqrt((1+tau**2)))
  
  return(
    list(var = Sigma,
         Y = Y,
         Z = Z, 
         f_Y = f_Y,
         g_Y = g_Y,
         pvals = pvals,
         pvals_fission = pvals_fission)
  )
}


summary.BH <- function(pvals, H0, vals, mu, var, alpha.list,
                       debug=FALSE) {
  n <- length(pvals)
  nfrej <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    sum(pvals[!H0] < alpha, na.rm = TRUE)
  })
  ntrej <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    sum(pvals[H0] < alpha, na.rm = TRUE)
  })
  true_mean <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    mean(mu[pvals < alpha], na.rm = TRUE)
  })
  measured_mean <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    mean(vals[pvals < alpha], na.rm = TRUE)
  })
  nrej <- nfrej + ntrej
  FDP <- nfrej / pmax(nrej, 1)
  power <- ntrej / max(sum(H0), 1)
  sd <- sqrt(var/nrej)
  df <- data.frame(nrej = nrej, 
                   FDP = FDP, 
                   power = power, 
                   mean=measured_mean, 
                   true_mean=true_mean,
                   sd=sd)
  return(df)
}


summary.AdaPT <- function(adapt, H0, pvals, vals, mu, var){
  nfrej <- apply(adapt$s, 2, function(s){
    bool.mask <- (pvals <= s)
    sum(bool.mask[!H0], na.rm = TRUE)
  })
  ntrej <- apply(adapt$s, 2, function(s){
    bool.mask <- (pvals <= s)        
    sum(bool.mask[H0], na.rm = TRUE)
  })
  true_mean <- apply(adapt$s, 2, function(s){
    bool.mask <- (pvals <= s)        
    mean(mu[bool.mask], na.rm = TRUE)
  })
  measured_mean <- apply(adapt$s, 2, function(s){
    bool.mask <- (pvals <= s)        
    mean(vals[bool.mask], na.rm = TRUE)
  })
  
  nrej <- nfrej + ntrej
  FDP <- nfrej / pmax(nrej, 1)
  power <- ntrej / max(sum(H0),1)
  sd <- sqrt(var/nrej)
  df <- data.frame(nrej = nrej, FDP = FDP, power = power, mean=measured_mean, true_mean =true_mean,sd=sd)
  return(df)
}


summary.STAR <- function(STAR.obj, H0, vals, mu, var){
  nfrej <- apply(STAR.obj$mask, 2, function(x){
    sum(x[!H0], na.rm = TRUE)
  })
  ntrej <- apply(STAR.obj$mask, 2, function(x){
    sum(x[H0], na.rm = TRUE)
  })
  true_mean <- apply(STAR.obj$mask, 2, function(x){
    mean(mu[x])
  })
  measured_mean <- apply(STAR.obj$mask, 2, function(x){
    mean(vals[x])
  })
  nrej <- nfrej + ntrej
  FDP <- nfrej / pmax(nrej, 1)
  power <- ntrej / max(sum(H0),1)
  sd <- sqrt(var/nrej)
  df <- data.frame(nrej = nrej, FDP = FDP, power = power, mean=measured_mean, true_mean =true_mean, sd=sd)
  return(df)
}


interactive_testing_exp <- function(idx, debug=FALSE) {
  if (debug) print(paste('Exp run', idx))
  
  # Generate the data
  out <- gen_data(tau)
  for (param in names(out)) {assign(param, out[[param]])}
  
  # Reject with BH, AdaPT, STAR for full dataset
  STAR.obj1 <- STAR.convex(pvals, x,
                           alpha.list = alpha.list,
                           type = "model-assist",
                           update.mask.params = list(dir="min"),
                           num.steps.update.score = num.steps.update.score,
                           score.params = score.params)
  AdaPT.obj1 <- AdaPT(x, pvals, cov.formula = cov.formula,
                      q.list = alpha.list, plot.quiet = TRUE)
  
  # Compute the results for the three CI procedures for the full dataset
  BH.result_full <- summary.BH(pvals, H0, Y, mu, 1, alpha.list, debug=debug)
  STAR.result_full <- summary.STAR(STAR.obj1, H0, Y, mu, 1)
  adapt.result_full <- summary.AdaPT(AdaPT.obj1, H0, pvals, Y, mu, 1)
  
  #Reject with BH, AdaPT, STAR for fissioned data
  STAR.obj2 <- STAR.convex(pvals_fission, x, 
                           alpha.list = alpha.list,
                           type = "model-assist",
                           update.mask.params = list(dir = "min"),
                           num.steps.update.score = num.steps.update.score,
                           score.params = score.params)
  AdaPT.obj2 <- AdaPT(x, pvals_fission, cov.formula = cov.formula,
                      q.list = alpha.list, plot.quiet = TRUE)
  
  # Compute the results for the three CI procedures for the fissioned dataset
  var_fission <- 1 + (1/(tau**2))
  BH.result.fission <- summary.BH(pvals_fission, H0, g_Y, mu, var_fission, alpha.list)
  STAR.result.fission <- summary.STAR(STAR.obj2, H0, g_Y, mu, var_fission)
  adapt.result.fission <- summary.AdaPT(AdaPT.obj2, H0, pvals, g_Y, mu, var_fission)
  
  return(list(x=x,mu=mu,tau=tau,H0=H0,mu=mu,alt=alt,null=null,var=var,var_fission=var_fission,
              Y=Y,f_Y=f_Y,g_Y=g_Y,pvals=pvals,pvals_fission=pvals_fission,
              BH.result.fission=BH.result.fission, BH.result.full=BH.result_full,
              adapt.result.fission=adapt.result.fission,adapt.result.full=adapt.result_full,
              STAR.result.fission=STAR.result.fission, STAR.result.full = STAR.result_full,
              STAR.full.object=STAR.obj1,AdaPT.full.object=AdaPT.obj1,
              STAR.mask.object=STAR.obj2,AdaPT.mask.object=AdaPT.obj2))
}


# Functions to generate plots
compute_trial_stats <- function(one.exp) {
  alpha.FCR <- 0.2
  alpha.list <- seq(0.01, 0.3, 0.02)
  
  q <- nrow(one.exp$x)
  pvals <- one.exp$pvals
  n <- length(pvals)
  pvals_fission <- one.exp$pvals_fission
  tau <- one.exp$tau
  
  BH_full_CI <- sapply(alpha.list, function(x) {
    khat <- max(c(0, which(sort(pvals) <= alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    x <- which(pvals < alpha)
    nrej <- length(x)
    if (nrej == 0) return(rep(NaN, 5))
    alpha.2 <- alpha.FCR*nrej/q
    ci_u <- one.exp$Y[x] + qnorm(1-alpha.FCR/2)
    ci_l <- one.exp$Y[x] - qnorm(1-alpha.FCR/2)
    cov <- (one.exp$mu[x] <= ci_u) & (one.exp$mu[x] >= ci_l)
    ret <- c(mean(one.exp$mu[x]),
             mean(one.exp$Y[x]),
             sum(cov),
             length(cov),
             mean(ci_u-ci_l))
    return(ret)
  })
  
  ada_full_CI <- apply(one.exp$AdaPT.full.object$s, 2, function(s) {
    return(1)
  })
  
  return(list(bh.f.ci = BH_full_CI))
}

gen_plot_df <- function(it_tau) {
  tmp <- (one.exp$pvals <= s)  
  x <- (one.exp$pvals_fission <= s)
}





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (tau in seq(0.1, 0.9, 0.1)) {
  load(paste0('abb_results/it_results/interactive_tau_',tau,'.Rdata'))
  assign(paste0('it_tau', tau*10), interactive_testing_res)
}




create_plots <- FALSE
if (create_plots) {

}
